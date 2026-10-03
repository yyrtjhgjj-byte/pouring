"""注ぎ方（戦略）を評価する。

戦略 = (コップの傾き θ, ノズルと縁の隙間 gap, 狙う場所 aim, レバーを開ける速さ τ_v)
評価指標:
  - 注ぎ始めの「水滴の柱」（Worthington ジェット）の高さ、ノズル到達の有無
  - 定常の注ぎで出る微細飛沫（気泡破裂ジェット液滴）の量と到達高さ
  - ノズルに飛沫が届く回数（衛生面: ノズルの汚染）
  - コップの縁がノズルに触れる確率（手ブレ）
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import math

import numpy as np
from scipy.stats import norm

from .props import Fluid, G
from .server import Server
from .cup import Cup, place_cup
from . import jet as jetmod
from . import startup
from . import bubbles
from . import droplets
from .wall import wall_film


@dataclass(frozen=True)
class Strategy:
    theta_deg: float = 0.0         # コップの傾き [deg]
    gap: float = 0.03              # ノズル先端〜縁の最上点の鉛直距離 [m]
    aim: str = "center"            # "center" or "wall"
    wall_offset: float = 0.01      # aim="wall" のとき、水面からどれだけ上の壁に当てるか [m]
    tau_valve: float | None = None  # レバー開栓時間 [s]（None ならサーバー既定）

    def label(self) -> str:
        a = "壁沿い" if self.aim == "wall" else "中央"
        return f"θ={self.theta_deg:.0f}° gap={self.gap*100:.1f}cm {a}"


@dataclass
class Scenario:
    server: Server = field(default_factory=Server)
    cup: Cup = field(default_factory=Cup)
    V_initial: float = 120e-6      # 最初から入っている水 [m^3]
    V_add: float = 150e-6          # 追加で注ぐ量 [m^3]
    r_nozzle: float = 0.007        # ノズル外半径 [m]
    tremor: float = 0.003          # 手ブレ（上下）の標準偏差 [m]
    r_burst: float = 0.012         # 気泡が水面で弾ける範囲の半径 [m]


@dataclass
class Evaluation:
    strategy: Strategy
    feasible: bool
    note: str
    mode: str                      # "plunge"（水面に直撃）/ "wall"（壁沿い）
    H_fall: float                  # 自由落下距離 [m]
    U_impact: float                # 衝突速度 [m/s]
    U_eff: float                   # 水面に突っ込む実効速度 [m/s]
    roughness: float               # 実効の表面乱れ
    V_onset: float                 # 巻き込み開始速度 [m/s]
    air_ratio: float               # 巻き込み空気/水
    drop_rate: float               # 飛沫発生率 [個/s]
    spray_h99: float               # 飛沫到達高さ 上位1% [m]
    H_surface_to_nozzle: float     # 水面からノズルまで [m]
    nozzle_hits: float             # 1 回の注ぎでノズルに当たる飛沫の期待個数
    rim_escape: float              # 1 回の注ぎで縁より上に飛ぶ飛沫の個数
    column_regime: str             # 注ぎ始めの柱のレジーム
    column_height: float           # 柱の到達高さ [m]（水面基準）
    column_hits_nozzle: bool
    head_Fr: float
    clearance: float               # 縁とノズルの最短距離 [m]
    p_contact: float               # 手ブレで縁がノズルに触れる確率
    badness: float = 0.0
    extras: dict = field(default_factory=dict)

    def as_row(self) -> dict:
        d = asdict(self)
        d.pop("extras")
        s = d.pop("strategy")
        d.update({f"s_{k}": v for k, v in s.items()})
        return d


def badness_score(ev: "Evaluation") -> float:
    """総合的な「ダメさ」（小さいほど良い）。重みは主観だが、すべての指標を見えるようにしておく。

    - ノズルに飛沫が当たる（衛生）: 1 個あたり 1 点
    - 注ぎ始めの柱がノズルに届く: 5 点、縁を越える: 1 点
    - 縁より上に飛ぶ微細飛沫: log10(1 + 個数)
    - 縁がノズルに触れる確率: ×20 点
    """
    b = ev.nozzle_hits
    b += 5.0 if ev.column_hits_nozzle else 0.0
    b += math.log10(1.0 + ev.rim_escape)
    b += 20.0 * ev.p_contact
    if not ev.feasible:
        b += 100.0
    return b


def evaluate(strategy: Strategy, sc: Scenario | None = None, V_now: float | None = None,
             with_head: bool = True) -> Evaluation:
    """ある水量 V_now（既定: 注ぎ始め）の瞬間について、戦略を評価する。"""
    sc = sc or Scenario()
    srv = sc.server if strategy.tau_valve is None else sc.server.with_(tau_valve=strategy.tau_valve)
    fl: Fluid = srv.fluid
    V = sc.V_initial if V_now is None else V_now
    theta = math.radians(strategy.theta_deg)
    pose, info = place_cup(sc.cup, theta, V, strategy.gap, aim=strategy.aim,
                           wall_offset=strategy.wall_offset)
    zw = info["z_water"]
    feasible, note = True, ""
    if info["spill"]:
        feasible, note = False, "この傾きでは水がこぼれる"

    # 水流の着地点
    mode = "plunge"
    Zp_hit = pose.left_wall_Zp_at_x(0.0)
    z_hit_wall = float(pose.left_wall_point(Zp_hit)[2]) if 0.0 <= Zp_hit <= sc.cup.height else -1e9
    if strategy.aim == "wall" and z_hit_wall > zw:
        mode = "wall"
        H_fall = -z_hit_wall
    else:
        H_fall = -zw
    if Zp_hit > sc.cup.height:
        feasible, note = False, "水流がコップの口に入らない"

    imp = jetmod.impact_state(srv, max(H_fall, 1e-4))
    U, a = imp["U"], imp["a"]
    extras = {"impact": imp, "pose": pose, "info": info}

    if mode == "wall":
        phi = pose.wall_angle
        s_path = (z_hit_wall - zw) / math.cos(phi)
        film = wall_film(U, a, srv.Q, phi, s_path, imp["roughness"], fl)
        U_eff, rough, d_eff = film.U_normal, film.roughness, max(film.h_entry * 2, 1e-4)
        burst_offset = film.offset_x
        extras["film"] = film
    else:
        U_eff, rough, d_eff = U, imp["roughness"], 2 * a
        burst_offset = 0.0

    # 定常の巻き込みと飛沫
    Ve = bubbles.onset_velocity(rough, fl)
    ar = bubbles.air_entrainment_ratio(U_eff, d_eff, H_fall, srv.d0, srv.v0, rough, fl)
    if mode == "wall":
        ar *= 0.3  # 裏側が壁で空洞ができない分、巻き込み効率を下げる（仮定）
    R_med = bubbles.bubble_size_median(U_eff, max(d_eff, 2 * a), fl)
    spray = bubbles.spray_from_bubbles(ar * srv.Q, R_med, fl)
    t_pour = sc.V_add / srv.Q
    H_s2n = -zw
    # ノズルに当たる確率（幾何）: 飛沫は破裂域にばらまかれ、ほぼ鉛直に上がる
    geo = (sc.r_nozzle / sc.r_burst) ** 2 * math.exp(-(burst_offset / sc.r_burst) ** 2)
    nozzle_hits = spray.rate_above(H_s2n) * min(geo, 1.0) * t_pour
    h_rim = info["lowest_rim"] - zw
    rim_escape = spray.rate_above(max(h_rim, 0.0)) * t_pour * 0.3

    # 注ぎ始め: 頭部の衝突 → Worthington ジェット（水滴の柱）
    col_regime, col_h, head_Fr, col_hit = "none", 0.0, 0.0, False
    if with_head:
        head = startup.simulate_head(srv, max(H_fall, 1e-3))
        head_Fr = head.froude()
        extras["head"] = head
        if mode == "plunge":
            Vj = startup.worthington_jet_velocity(head.U, head_Fr)
            col_regime = startup.worthington_regime(head_Fr)
            ratio = Vj / max(head.U, 1e-6)
            r_col = head.D * (0.25 - 0.2 * float(np.clip((ratio - 0.7) / 3.05, 0, 1)))
            col_h = float(droplets.max_height(r_col, Vj, rho_l=fl.rho)) if Vj > 0.05 else Vj ** 2 / (2 * G)
            col_hit = col_h >= H_s2n * 0.98
            extras["column"] = {"Vj": Vj, "r_drop": r_col}
        else:
            col_regime = "wall"   # 頭部は壁に当たって膜に均される → 柱は立たない

    clearance = pose.nozzle_clearance(sc.r_nozzle)
    p_contact = float(norm.cdf(-clearance / sc.tremor))

    ev = Evaluation(
        strategy=strategy, feasible=feasible, note=note, mode=mode, H_fall=H_fall,
        U_impact=U, U_eff=U_eff, roughness=rough, V_onset=Ve, air_ratio=ar,
        drop_rate=spray.drop_rate, spray_h99=spray.h_drop_p99,
        H_surface_to_nozzle=H_s2n, nozzle_hits=nozzle_hits, rim_escape=rim_escape,
        column_regime=col_regime, column_height=col_h, column_hits_nozzle=col_hit,
        head_Fr=head_Fr, clearance=clearance, p_contact=p_contact, extras=extras,
    )
    extras["spray"] = spray
    ev.badness = badness_score(ev)
    return ev


def evaluate_pour(strategy: Strategy, sc: Scenario | None = None, n_steps: int = 6) -> dict:
    """注ぎ始めから注ぎ終わりまでの複数時点で評価し、まとめる（傾き・隙間は一定のまま）。"""
    sc = sc or Scenario()
    Vs = np.linspace(sc.V_initial, sc.V_initial + sc.V_add, n_steps)
    evs = [evaluate(strategy, sc, V_now=V, with_head=(i == 0)) for i, V in enumerate(Vs)]
    first = evs[0]
    frac = 1.0 / n_steps
    return {
        "strategy": strategy,
        "feasible": all(e.feasible for e in evs),
        "nozzle_hits": sum(e.nozzle_hits for e in evs) * frac,
        "rim_escape": sum(e.rim_escape for e in evs) * frac,
        "column_hits_nozzle": first.column_hits_nozzle,
        "column_height": first.column_height,
        "p_contact": max(e.p_contact for e in evs),
        "evals": evs,
    }
