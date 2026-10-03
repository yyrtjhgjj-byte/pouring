"""半経験モデルのパラメータ掃引 → results/model.json（レポートのグラフ用）。

  python3 scripts/sweep_model.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pourphys import Scenario, Strategy, evaluate, Server, COLD, ROOM, HOT  # noqa: E402
from pourphys import jet as jetmod, bubbles, startup, droplets  # noqa: E402
from pourphys.optimize import tilt_schedule  # noqa: E402
from pourphys.wall import wall_film  # noqa: E402
from pourphys.props import G  # noqa: E402


def plunge_curve(srv: Server, Hs, r_nozzle=0.007, r_burst=0.012, t_pour=6.0):
    """水面に真っすぐ落とすときの物理量を、落下距離 H の関数として。"""
    fl = srv.fluid
    rows = []
    for H in Hs:
        imp = jetmod.impact_state(srv, H)
        U, a, rough = imp["U"], imp["a"], imp["roughness"]
        Ve = bubbles.onset_velocity(rough, fl)
        ar = bubbles.air_entrainment_ratio(U, 2 * a, H, srv.d0, srv.v0, rough, fl)
        Rm = bubbles.bubble_size_median(U, 2 * a, fl)
        sp = bubbles.spray_from_bubbles(ar * srv.Q, Rm, fl)
        geo = (r_nozzle / r_burst) ** 2
        rows.append({"H": H, "U": U, "Ve": Ve, "rough": rough, "air": ar, "R_med": Rm,
                     "bubble_rate": sp.bubble_rate, "drop_rate": sp.drop_rate,
                     "h50": sp.h_drop_median, "h99": sp.h_drop_p99,
                     "nozzle_hits_per_pour": sp.rate_above(H) * geo * t_pour,
                     "d_jet": 2 * a})
    return rows


def wall_curve(srv: Server, Hs, phi_deg: float, offset=0.01):
    """内壁（鉛直から phi 傾いた面）に当てて、offset 下の水面に膜で入る場合。H は壁に当たるまで。"""
    fl = srv.fluid
    phi = math.radians(phi_deg)
    rows = []
    for H in Hs:
        imp = jetmod.impact_state(srv, H)
        s_path = offset / math.cos(phi)
        film = wall_film(imp["U"], imp["a"], srv.Q, phi, s_path, imp["roughness"], fl)
        Ve = bubbles.onset_velocity(film.roughness, fl)
        ar = 0.3 * bubbles.air_entrainment_ratio(film.U_normal, 2 * film.h_entry, H, srv.d0, srv.v0, film.roughness, fl)
        Rm = bubbles.bubble_size_median(film.U_normal, 2 * imp["a"], fl)
        sp = bubbles.spray_from_bubbles(ar * srv.Q, Rm, fl)
        rows.append({"H": H, "U": imp["U"], "Ueff": film.U_normal, "Ve": Ve, "air": ar,
                     "bubble_rate": sp.bubble_rate, "drop_rate": sp.drop_rate, "splash_K": film.splash_K})
    return rows


def crossing(rows, key_u="U", key_v="Ve"):
    for a, b in zip(rows[:-1], rows[1:]):
        da, db = a[key_u] - a[key_v], b[key_u] - b[key_v]
        if da <= 0 < db:
            return a["H"] + (b["H"] - a["H"]) * (-da) / (db - da)
    return None


def main():
    out = {}
    srv = Server()
    Hs = np.round(np.arange(0.02, 0.301, 0.005), 4)
    out["server"] = {"Q_Lmin": srv.Q * 60000, "d0_mm": srv.d0 * 1e3, "v0": srv.v0, "noise": srv.noise}
    out["plunge"] = plunge_curve(srv, Hs)
    out["H_star"] = crossing(out["plunge"])
    out["wall"] = {str(p): wall_curve(srv, Hs, p) for p in (4, 15, 30, 40)}
    # 水温
    out["temps"] = {}
    for name, fl in (("cold", COLD), ("room", ROOM), ("hot", HOT)):
        s2 = srv.with_(fluid=fl)
        rows = plunge_curve(s2, Hs)
        out["temps"][name] = {"T": fl.T_c, "sigma": fl.sigma, "mu": fl.mu, "H_star": crossing(rows),
                              "rows": [{"H": r["H"], "U": r["U"], "Ve": r["Ve"], "air": r["air"]} for r in rows]}
    # ノズルでの乱れ（サーバーの個体差・ボトルの「ゴボッ」）
    out["noise"] = []
    for eps in (0.005, 0.01, 0.02, 0.04, 0.08):
        rows = plunge_curve(srv.with_(noise=eps), Hs)
        out["noise"].append({"noise": eps, "H_star": crossing(rows)})
    # 流量
    out["flow"] = []
    for q in (1.0, 1.5, 2.0, 2.5):
        s2 = srv.with_(Q=q / 60000)
        rows = plunge_curve(s2, Hs)
        out["flow"].append({"Q_Lmin": q, "v0": s2.v0, "H_star": crossing(rows)})
    # レバーの開け方 → 先頭の水塊
    out["lever"] = []
    for tau in (0.02, 0.04, 0.06, 0.1, 0.15, 0.25, 0.4):
        for H in (0.06, 0.10, 0.15):
            h = startup.simulate_head(srv.with_(tau_valve=tau), H)
            out["lever"].append({"tau": tau, "H": H, "mass_g": h.mass * 1e3, "U": h.U, "D_mm": h.D * 1e3,
                                 "Fr": h.froude(), "kinematic": h.kinematic})
    # 水量ごとの傾き
    sc = Scenario()
    sched = tilt_schedule(sc, np.arange(40e-6, 300e-6, 10e-6))
    out["schedule"] = [{k: (v * 1e6 if k == "V" else v) for k, v in r.items()} for r in sched]
    out["cup"] = {"volume_mL": sc.cup.volume * 1e6, "r_bottom_mm": sc.cup.r_bottom * 1e3,
                  "r_top_mm": sc.cup.r_top * 1e3, "height_mm": sc.cup.height * 1e3}
    # 戦略の比較（120 mL から 150 mL 足す）
    comps = {
        "遠くから中央（隙間12cm）": Strategy(0, 0.12, "center"),
        "ふつう（隙間5cm・中央）": Strategy(0, 0.05, "center"),
        "近づけて中央（隙間1cm）": Strategy(0, 0.01, "center"),
        "直立のまま内壁をかすめる（隙間3cm）": Strategy(0, 0.03, "wall"),
        "30°傾けて内壁（隙間2.5cm）": Strategy(30, 0.025, "wall"),
        "30°傾けて内壁・そっと開栓": Strategy(30, 0.025, "wall", tau_valve=0.25),
    }
    out["strategies"] = []
    for name, st in comps.items():
        ev = evaluate(st, sc)
        out["strategies"].append({
            "name": name, "mode": ev.mode, "H": ev.H_fall, "U": ev.U_impact, "Ueff": ev.U_eff,
            "Ve": ev.V_onset, "air": ev.air_ratio, "drop_rate": ev.drop_rate,
            "nozzle_hits": ev.nozzle_hits, "rim_escape": ev.rim_escape, "clearance_mm": ev.clearance * 1e3,
            "p_contact": ev.p_contact, "head_Fr": ev.head_Fr, "column": ev.column_regime,
            "column_h": ev.column_height, "badness": ev.badness})
    # 気泡 → ジェット液滴（説明用）
    R = np.geomspace(20e-6, 1.6e-3, 40)
    rd, Vd, n = bubbles.jet_drops(R, ROOM)
    hh = droplets.max_height(rd, Vd)
    out["jetdrops"] = [{"R_um": float(a * 1e6), "rd_um": float(b * 1e6), "Vd": float(c), "n": int(d),
                        "h_cm": float(e * 100), "h_vac_cm": float(c * c / (2 * G) * 100)}
                       for a, b, c, d, e in zip(R, rd, Vd, n, hh)]
    (ROOT / "results").mkdir(exist_ok=True)
    with open(ROOT / "results" / "model.json", "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    print("H* (room) =", out["H_star"])
    for name, t in out["temps"].items():
        print(name, "H* =", t["H_star"])
    for r in out["noise"]:
        print("noise", r)
    for r in out["flow"]:
        print("flow", r)
    for s in out["strategies"]:
        print(s["name"], {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items() if k != "name"})


if __name__ == "__main__":
    main()
