"""注ぎ方の探索：どの傾き・隙間・狙いが一番マシか。

- grid_search: (傾き, 隙間, 狙い, 開栓時間) の格子を総当たりで評価
- pareto: 「ノズル汚染」と「飛沫の量」のトレードオフの非劣解
- tilt_schedule: コップが満ちていくにつれて、傾きをどう戻すべきか
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from .cup import place_cup
from .strategy import Strategy, Scenario, evaluate


def grid_search(sc: Scenario, thetas, gaps, aims=("center", "wall"), taus=(None,), V=None,
                with_head=True):
    out = []
    for th, gap, aim, tau in itertools.product(thetas, gaps, aims, taus):
        if aim == "center" and th > 0:
            continue
        ev = evaluate(Strategy(th, gap, aim, tau_valve=tau), sc, V_now=V, with_head=with_head)
        out.append(ev)
    return out


def pareto(evals, keys=("nozzle_hits", "drop_rate")):
    """keys の全てで劣らない評価だけを残す（小さいほど良い）。"""
    pts = np.array([[getattr(e, k) for k in keys] for e in evals])
    keep = []
    for i, p in enumerate(pts):
        dominated = np.any(np.all(pts <= p, axis=1) & np.any(pts < p, axis=1))
        if not dominated:
            keep.append(evals[i])
    return keep


def wall_aim_feasible(sc: Scenario, theta_deg: float, V: float, gap: float,
                      offset: float = 0.01, rim_margin: float = 0.015, spill_margin: float = 0.008):
    """壁沿いに注げるか: こぼれない & 当たる点が水面より上で縁より十分下。"""
    th = math.radians(theta_deg)
    pose, info = place_cup(sc.cup, th, V, gap, aim="wall", wall_offset=offset)
    zw = info["z_water"]
    if zw > info["lowest_rim"] - spill_margin:
        return False, pose, info
    Zp_hit = pose.left_wall_Zp_at_x(0.0)
    if not (0 < Zp_hit < sc.cup.height - rim_margin / max(math.cos(pose.wall_angle), 0.3)):
        return False, pose, info
    z_hit = float(pose.left_wall_point(Zp_hit)[2])
    if z_hit < zw + 0.003:
        return False, pose, info
    return True, pose, info


def tilt_schedule(sc: Scenario, volumes, thetas=np.arange(0, 51, 1), gap=0.025, offset=0.01):
    """水量ごとに、壁沿いに注げる最大の傾きと、そのときの評価。"""
    rows = []
    for V in volumes:
        best = None
        for th in thetas[::-1]:
            ok, pose, info = wall_aim_feasible(sc, th, V, gap, offset)
            if ok:
                best = th
                break
        if best is None:
            rows.append({"V": V, "theta_max": None})
            continue
        ev = evaluate(Strategy(best, gap, "wall", wall_offset=offset), sc, V_now=V, with_head=False)
        # 推奨: 最大傾きから 5° 余裕をみる（手ブレ・揺れ）
        rec = max(best - 5, 0)
        ev_rec = evaluate(Strategy(rec, gap, "wall", wall_offset=offset), sc, V_now=V, with_head=False)
        rows.append({"V": V, "theta_max": float(best), "theta_rec": float(rec),
                     "Ueff_max": ev.U_eff, "Ueff_rec": ev_rec.U_eff, "Ve": ev_rec.V_onset,
                     "mode_rec": ev_rec.mode, "air_rec": ev_rec.air_ratio, "H_rec": ev_rec.H_fall})
    return rows
