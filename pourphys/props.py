"""水と空気の物性値（温度依存）。

ウォーターサーバーには冷水（~6 °C）と温水（~85 °C）があるので、
粘度・表面張力・密度を温度の関数として持っておく。
温水は粘度も表面張力も小さいので、同じ注ぎ方でも飛沫が出やすい。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

G = 9.80665  # 重力加速度 [m/s^2]


def water_density(T_c: float) -> float:
    """水の密度 [kg/m^3]（Tanaka et al. 2001 型の経験式, 0–100 °C で誤差 ~0.1%）。"""
    a1, a2, a3, a4, a5 = -3.983035, 301.797, 522528.9, 69.34881, 999.974950
    return a5 * (1.0 - (T_c + a1) ** 2 * (T_c + a2) / (a3 * (T_c + a4)))


def water_viscosity(T_c: float) -> float:
    """水の粘度 [Pa s]（Vogel 式）。"""
    T = T_c + 273.15
    return 2.414e-5 * 10.0 ** (247.8 / (T - 140.0))


def water_surface_tension(T_c: float) -> float:
    """水の表面張力 [N/m]（IAPWS 1994）。"""
    Tc = 647.096
    tau = 1.0 - (T_c + 273.15) / Tc
    return 235.8e-3 * tau ** 1.256 * (1.0 - 0.625 * tau)


def air_viscosity(T_c: float) -> float:
    """空気の粘度 [Pa s]（Sutherland 式）。"""
    T = T_c + 273.15
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


@dataclass(frozen=True)
class Fluid:
    """液体（水）と周囲気体（空気）の物性セット。"""

    T_c: float = 20.0
    rho: float = 998.2
    mu: float = 1.0016e-3
    sigma: float = 0.07275
    rho_air: float = 1.204
    mu_air: float = 1.81e-5

    @classmethod
    def water(cls, T_c: float = 20.0) -> "Fluid":
        # 空気側は室温 (25 °C) 固定。湯気の効果などは無視。
        return cls(
            T_c=T_c,
            rho=water_density(T_c),
            mu=water_viscosity(T_c),
            sigma=water_surface_tension(T_c),
            rho_air=1.184,
            mu_air=air_viscosity(25.0),
        )

    # --- 派生量 ---
    @property
    def nu(self) -> float:
        return self.mu / self.rho

    @property
    def capillary_length(self) -> float:
        """毛管長 l_c = sqrt(σ/(ρg)) [m]。水で ~2.7 mm。"""
        return math.sqrt(self.sigma / (self.rho * G))

    @property
    def visco_capillary_length(self) -> float:
        """粘性-毛管長 l_μ = μ²/(ρσ) [m]。水で ~14 nm。"""
        return self.mu ** 2 / (self.rho * self.sigma)

    @property
    def visco_capillary_velocity(self) -> float:
        """V_μ = σ/μ [m/s]。水で ~73 m/s。"""
        return self.sigma / self.mu

    def weber(self, U: float, L: float) -> float:
        return self.rho * U * U * L / self.sigma

    def reynolds(self, U: float, L: float) -> float:
        return self.rho * U * L / self.mu

    def froude(self, U: float, L: float) -> float:
        """Fr = U²/(gL)（落下液滴の文献流儀。平方根を取らない方）。"""
        return U * U / (G * L)

    def ohnesorge(self, L: float) -> float:
        return self.mu / math.sqrt(self.rho * self.sigma * L)

    def bond(self, L: float) -> float:
        return self.rho * G * L * L / self.sigma


COLD = Fluid.water(6.0)    # 冷水
ROOM = Fluid.water(20.0)   # 常温水
HOT = Fluid.water(85.0)    # 温水

PRESETS = {"cold": COLD, "room": ROOM, "hot": HOT}
