"""ウォーターサーバー（吐出側）のモデル。

- 定常吐出: 流量 Q と吐出口径 d0 から初速 v0 を決める。
- 開栓過渡: レバーを開く速さ（開栓時間 τ_v）で流量が立ち上がる。
- 脈動: ボトル式サーバーの「ゴボッ」（空気がボトルに入る）による流量の周期変動。

数値は国内の一般的なサーバー仕様（冷水 1.5–2 L/min 程度）を想定した既定値。
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

import numpy as np

from .props import Fluid, ROOM


def smoothstep(s):
    s = np.clip(s, 0.0, 1.0)
    return s * s * (3.0 - 2.0 * s)


@dataclass(frozen=True)
class Server:
    """吐出口（ノズル）の条件。"""

    Q: float = 1.5e-3 / 60.0     # 定常流量 [m^3/s]（1.5 L/min）
    d0: float = 7.0e-3           # 吐出口の水柱直径 [m]
    tau_valve: float = 0.06      # レバー開栓にかかる時間 [s]（ガッと開けると ~0.05 s）
    pulsation: float = 0.03      # 流量脈動の相対振幅（ボトル式の「ゴボッ」）
    f_pulse: float = 3.0         # 脈動周波数 [Hz]
    noise: float = 0.02          # 吐出口での水柱表面の初期乱れ（半径比）
    fluid: Fluid = field(default=ROOM)

    @property
    def area(self) -> float:
        return math.pi * self.d0 ** 2 / 4.0

    @property
    def v0(self) -> float:
        """定常時の吐出速度 [m/s]。"""
        return self.Q / self.area

    def flow(self, t):
        """時刻 t [s] での流量 [m^3/s]（開栓立ち上がり＋脈動）。"""
        t = np.asarray(t, dtype=float)
        ramp = smoothstep(t / self.tau_valve) if self.tau_valve > 0 else (t >= 0).astype(float)
        puls = 1.0 + self.pulsation * np.sin(2.0 * math.pi * self.f_pulse * t)
        return self.Q * ramp * puls

    def exit_velocity(self, t):
        """時刻 t での吐出速度 [m/s]。口径一杯に流れる（満管）と仮定。"""
        return self.flow(t) / self.area

    def with_(self, **kw) -> "Server":
        return replace(self, **kw)


# 代表的なプリセット
PRESET_SERVERS = {
    "standard": Server(),
    "fast_flow": Server(Q=2.5e-3 / 60.0, d0=8.0e-3),
    "slow_flow": Server(Q=1.0e-3 / 60.0, d0=6.0e-3),
}
