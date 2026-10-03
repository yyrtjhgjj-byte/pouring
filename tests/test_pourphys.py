"""pourphys の基本的な性質のテスト（物理的に当然成り立つべきこと）。

  python3 -m pytest -q tests
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pourphys import props, jet, bubbles, droplets, startup  # noqa: E402
from pourphys import Server, Cup, Pose, place_cup, Strategy, Scenario, evaluate  # noqa: E402
from pourphys.optimize import tilt_schedule  # noqa: E402


def test_water_properties_match_tables():
    f = props.Fluid.water(20.0)
    assert f.rho == pytest.approx(998.2, abs=0.2)
    assert f.mu == pytest.approx(1.002e-3, rel=0.02)
    assert f.sigma == pytest.approx(0.0728, rel=0.01)
    hot = props.Fluid.water(85.0)
    assert hot.mu < f.mu / 2 and hot.sigma < f.sigma
    assert f.capillary_length == pytest.approx(2.7e-3, rel=0.03)


def test_falling_jet_conserves_mass_and_accelerates():
    srv = Server()
    js = jet.falling_jet(srv, z_max=0.3)
    Q = math.pi * js.a ** 2 * js.v
    assert np.allclose(Q, srv.Q, rtol=1e-6)
    assert np.all(np.diff(js.v) > 0)
    # 表面張力の補正は小さい: 自由落下の速さとの差は 1% 以内
    v_ff = np.sqrt(srv.v0 ** 2 + 2 * props.G * js.z)
    assert np.max(np.abs(js.v / v_ff - 1)) < 0.03
    assert np.all(np.diff(js.gain) >= -1e-12)


def test_rayleigh_plateau_most_unstable_wavenumber():
    f = props.ROOM
    a = 2e-3
    x = np.linspace(0.05, 0.99, 400)
    w = jet.rp_growth_rate(x, a, f)
    # 非粘性の最大成長は x ≈ 0.697
    assert x[np.argmax(w)] == pytest.approx(0.697, abs=0.02)
    assert np.all(jet.rp_growth_rate(np.array([1.0, 1.5]), a, f) <= 1e-12)


def test_cup_volume_and_tilted_water_level():
    cup = Cup()
    assert cup.volume * 1e6 == pytest.approx(335, abs=3)
    pose = Pose(cup, 0.0, 0.0, 0.0)
    V = 120e-6
    zw = pose.water_level(V)
    assert zw == pytest.approx(cup.upright_level(V), abs=2e-4)
    # 傾けても体積は変わらない
    p2 = Pose(cup, math.radians(30), 0.0, 0.0)
    zw2 = p2.water_level(V)
    assert p2.water_volume_below(zw2) == pytest.approx(V, rel=1e-3)
    # 傾けるとこぼれずに入る量は減る
    assert p2.max_volume() < cup.volume


def test_jet_drop_scaling_reasonable():
    f = props.ROOM
    R = np.array([0.1e-3, 0.25e-3, 0.5e-3])
    rd, Vd, n = bubbles.jet_drops(R, f)
    # 「1/10 則」: 先頭液滴は気泡直径のおよそ 1/10
    assert np.all((rd / R > 0.05) & (rd / R < 0.2))
    # 小さい泡ほど速い
    assert np.all(np.diff(Vd) < 0)
    assert np.all(n > 0)


def test_drag_limits_small_drop_height():
    h = droplets.max_height(20e-6, 5.0)
    assert h < 0.02            # 20 µm の粒は 5 m/s でも 2 cm 届かない
    assert droplets.max_height(1e-3, 1.0) == pytest.approx(1.0 / (2 * props.G), rel=0.1)


def test_entrainment_increases_with_fall_height():
    srv = Server()
    ars = []
    for H in (0.04, 0.08, 0.12, 0.2):
        imp = jet.impact_state(srv, H)
        ars.append(bubbles.air_entrainment_ratio(imp["U"], imp["d"], H, srv.d0, srv.v0, imp["roughness"], srv.fluid))
    assert all(b >= a for a, b in zip(ars, ars[1:]))
    assert ars[-1] > 100 * max(ars[0], 1e-9)


def test_tilting_and_wall_reduce_spray():
    sc = Scenario()
    straight = evaluate(Strategy(0, 0.05, "center"), sc, with_head=False)
    tilted = evaluate(Strategy(30, 0.025, "wall"), sc, with_head=False)
    assert tilted.mode == "wall"
    assert tilted.U_eff < straight.U_eff
    assert tilted.drop_rate < straight.drop_rate / 100


def test_tilt_schedule_decreases_with_volume():
    sc = Scenario()
    rows = [r for r in tilt_schedule(sc, np.arange(60e-6, 220e-6, 40e-6)) if r["theta_max"] is not None]
    th = [r["theta_max"] for r in rows]
    assert th == sorted(th, reverse=True)


def test_startup_head_forms_only_for_fast_valve():
    fast = startup.simulate_head(Server(tau_valve=0.02), 0.15)
    slow = startup.simulate_head(Server(tau_valve=0.4), 0.15)
    assert fast.kinematic and not slow.kinematic
