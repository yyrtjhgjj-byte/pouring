"""ノズルから落ちる水柱（自由落下ジェット）のモデル。

- 重力で加速しながら細くなる（質量保存）。表面張力によるわずかな圧力も含めた Bernoulli 式。
- Rayleigh–Plateau 不安定: 水柱表面の小さなくびれが、落下中に指数的に成長する。
  水柱は加速で引き伸ばされるので、波長も伸び、半径も縮む。それを物質要素に沿って積分する。
  → 着水時の「表面の凸凹（乱れ）」と、水柱が粒に千切れる距離（分裂長）が出る。

z はノズルから下向きに測った落下距離 [m]。
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from functools import lru_cache

import numpy as np
from scipy.special import i0e, i1e

from .props import G
from .server import Server


@dataclass
class JetState:
    z: np.ndarray        # 落下距離 [m]
    v: np.ndarray        # 速度 [m/s]
    a: np.ndarray        # 半径 [m]
    t: np.ndarray        # ノズルからの経過時間（物質要素）[s]
    gain: np.ndarray     # 最も危険な擾乱の成長率積分 ln(A/A0)
    k0_star: np.ndarray  # その擾乱のノズルでの波数 [1/m]

    def at(self, z: float) -> dict:
        """落下距離 z での状態を補間して返す。"""
        f = lambda arr: float(np.interp(z, self.z, arr))
        return {"z": z, "v": f(self.v), "a": f(self.a), "t": f(self.t),
                "gain": f(self.gain), "k0_star": f(self.k0_star)}


def bernoulli_profile(server: Server, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """定常落下ジェットの速度 v(z) と半径 a(z)。

    1/2 v^2 + σ/(ρ a) = 1/2 v0^2 + σ/(ρ a0) + g z,   Q = π a^2 v
    を固定点反復で解く（表面張力項は低速時に数%効く程度）。
    """
    fl = server.fluid
    Q, a0, v0 = server.Q, server.d0 / 2.0, server.v0
    E = 0.5 * v0 ** 2 + fl.sigma / (fl.rho * a0) + G * z
    v = np.sqrt(v0 ** 2 + 2 * G * z)
    for _ in range(50):
        a = np.sqrt(Q / (math.pi * v))
        v_new = np.sqrt(np.maximum(2.0 * (E - fl.sigma / (fl.rho * a)), 1e-12))
        if np.max(np.abs(v_new - v)) < 1e-12:
            break
        v = v_new
    a = np.sqrt(Q / (math.pi * v))
    return v, a


def rp_growth_rate(x, a, fluid) -> np.ndarray:
    """Rayleigh–Plateau 成長率 ω [1/s]。

    非粘性 Rayleigh の分散関係 ω_R^2 = σ/(ρ a^3) x (1-x^2) I1(x)/I0(x) に、
    粘性による減速を Chandrasekhar 型の近似 ω = sqrt(ω_R^2 + s^2) - s,
    s = 3 ν x^2 / (2 a^2) で入れる。x = k a ≥ 1 は安定（成長 0 とする）。
    """
    x = np.asarray(x, dtype=float)
    cap = fluid.sigma / (fluid.rho * a ** 3)
    xs = np.clip(x, 1e-9, None)
    ratio = i1e(xs) / i0e(xs)  # 指数スケール版で桁あふれしない
    wR2 = np.where(x < 1.0, cap * xs * (1.0 - xs ** 2) * ratio, 0.0)
    s = 1.5 * fluid.nu * xs ** 2 / a ** 2
    return np.sqrt(wR2 + s ** 2) - s


def falling_jet(server: Server, z_max: float = 0.4, n: int = 800, n_k: int = 160) -> JetState:
    """ノズルから z_max まで、加速・細り・RP 成長を計算する。

    物質要素に沿った成長: ln(A/A0) = ∫ ω(k(t) a(t)) dt
      k(t) = k0 v0 / v(t)  （引き伸ばしで波長は v に比例して伸びる）
    ノズルでの波数 k0 を走査し、各 z で最大の成長を与える k0 を「最も危険な擾乱」とする。
    """
    fl = server.fluid
    z = np.linspace(0.0, z_max, n)
    v, a = bernoulli_profile(server, z)
    # 経過時間 t(z) = ∫ dz / v
    dz = np.diff(z)
    t = np.concatenate([[0.0], np.cumsum(dz * 0.5 * (1.0 / v[1:] + 1.0 / v[:-1]))])
    a0, v0 = a[0], v[0]
    # 候補波数（ノズル基準で x0 = k0 a0 ∈ [0.05, 1.5]）
    x0 = np.linspace(0.05, 1.5, n_k)
    k0 = x0 / a0
    # 各 (z, k0) で x = k a
    kz = k0[None, :] * (v0 / v[:, None])
    x = kz * a[:, None]
    omega = rp_growth_rate(x, a[:, None], fl)
    dt = np.diff(t)[:, None]
    gains = np.vstack([np.zeros((1, n_k)), np.cumsum(0.5 * (omega[1:] + omega[:-1]) * dt, axis=0)])
    idx = np.argmax(gains, axis=1)
    gain = gains[np.arange(n), idx]
    return JetState(z=z, v=v, a=a, t=t, gain=gain, k0_star=k0[idx])


@lru_cache(maxsize=64)
def _jet_cached(server: Server) -> JetState:
    return falling_jet(server, z_max=0.6)


def impact_state(server: Server, H: float) -> dict:
    """落下距離 H で着水する時の水柱の状態。

    返り値:
      U: 着水速度, a: 半径, d: 直径, roughness: 表面乱れ（半径比 δ/a）,
      broken: 着水前に粒に分裂しているか, breakup_length: 分裂長
    """
    js = _jet_cached(server) if H <= 0.6 else falling_jet(server, z_max=H * 1.5)
    st = js.at(H)
    # 乱れの振幅は ノズルでの初期乱れ × exp(gain)
    rough = min(server.noise * math.exp(st["gain"]), 1.0)
    # 分裂長: 乱れが半径に達する距離
    amp = server.noise * np.exp(js.gain)
    i_br = np.argmax(amp >= 1.0) if np.any(amp >= 1.0) else None
    L_br = float(js.z[i_br]) if i_br is not None else float("inf")
    return {
        "H": H, "U": st["v"], "a": st["a"], "d": 2 * st["a"],
        "t_fall": st["t"], "gain": st["gain"], "roughness": rough,
        "broken": H >= L_br, "breakup_length": L_br,
    }
