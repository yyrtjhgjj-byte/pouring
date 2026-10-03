"""飛沫（小さな水滴）の弾道計算：空気抵抗つき。

飛沫の到達高さは「打ち出し速度」だけでなく「大きさ」で決まる。
100 μm 以下の粒は空気抵抗ですぐ止まるので、速くても数 cm しか上がらない。

抵抗係数は Schiller–Naumann: C_D = 24/Re (1 + 0.15 Re^0.687)（Re < 1000）。
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

from .props import G

RHO_AIR = 1.184
MU_AIR = 1.85e-5


def _drag_k(vmag, r, rho_l):
    Re = 2.0 * r * vmag * RHO_AIR / MU_AIR
    Cd = 24.0 / max(Re, 1e-12) * (1.0 + 0.15 * Re ** 0.687) if Re < 1000.0 else 0.44
    # F = 1/2 ρa Cd A v^2,  m = ρl 4/3 π r^3  ->  |a| = 3 ρa Cd v^2 / (8 ρl r)
    return 3.0 * RHO_AIR * Cd * vmag / (8.0 * rho_l * r)


def fly(r: float, V: float, angle_deg: float = 90.0, rho_l: float = 998.0,
        z_floor: float = -1.0, t_max: float = 3.0):
    """半径 r [m], 速さ V [m/s], 仰角 angle_deg で打ち出した粒の軌道。

    返り値 dict: zmax（最高到達高さ）, x_at_zmax, t_apex, sol（solve_ivp の解）
    """
    th = math.radians(angle_deg)

    def rhs(t, y):
        x, z, vx, vz = y
        vmag = math.hypot(vx, vz) + 1e-15
        k = _drag_k(vmag, r, rho_l)
        return [vx, vz, -k * vx, -k * vz - G]

    def apex(t, y):
        return y[3]
    apex.terminal = False
    apex.direction = -1

    def floor(t, y):
        return y[1] - z_floor
    floor.terminal = True
    floor.direction = -1

    sol = solve_ivp(rhs, (0.0, t_max), [0.0, 0.0, V * math.cos(th), V * math.sin(th)],
                    method="LSODA", events=(apex, floor), rtol=1e-7, atol=1e-10,
                    dense_output=True)
    if sol.t_events[0].size:
        ta = sol.t_events[0][0]
        ya = sol.sol(ta)
        zmax, xa = float(ya[1]), float(ya[0])
    else:
        ta, zmax, xa = 0.0, 0.0, 0.0
    return {"zmax": zmax, "x_at_zmax": xa, "t_apex": ta, "sol": sol}


_CACHE_DIR = Path(__file__).resolve().parent / ".cache"


@lru_cache(maxsize=4)
def height_table(rho_l: float = 998.0):
    """最高到達高さ h(r, V) の補間表（鉛直打ち上げ）。log-log 補間。初回のみ計算してキャッシュ。"""
    rs = np.geomspace(2e-6, 4e-3, 40)
    Vs = np.geomspace(0.05, 40.0, 40)
    cache = _CACHE_DIR / f"height_table_{int(rho_l)}.npy"
    if cache.exists():
        H = np.load(cache)
    else:
        H = np.zeros((len(rs), len(Vs)))
        for i, r in enumerate(rs):
            for j, V in enumerate(Vs):
                H[i, j] = fly(r, V, rho_l=rho_l)["zmax"]
        _CACHE_DIR.mkdir(exist_ok=True)
        np.save(cache, H)
    return RegularGridInterpolator((np.log(rs), np.log(Vs)), np.log(np.maximum(H, 1e-12)),
                                   bounds_error=False, fill_value=None)


def max_height(r, V, rho_l: float = 998.0):
    """鉛直に打ち上げた粒の最高到達高さ [m]（補間表）。配列可。"""
    r = np.asarray(r, float)
    V = np.asarray(V, float)
    interp = height_table(float(round(rho_l)))
    lr, lV = np.broadcast_arrays(np.log(np.clip(r, 2e-6, 4e-3)), np.log(np.clip(V, 0.05, 40.0)))
    h = np.exp(interp(np.stack([lr, lV], axis=-1))).reshape(lr.shape)
    out = np.where(V <= 0.05, V ** 2 / (2 * G), h)
    return float(out) if out.ndim == 0 else out


def vacuum_height(V):
    """空気抵抗なしの到達高さ V²/2g。"""
    return np.asarray(V) ** 2 / (2 * G)
