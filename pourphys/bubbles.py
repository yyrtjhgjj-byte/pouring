"""空気の巻き込み → 気泡 → 水面で破裂 → 微細飛沫（ジェット液滴）。

「水面で微細な飛沫が上がる」の正体はほぼこれ。
着水した水柱が空気を巻き込むと、0.1–1 mm 程度の気泡が大量にでき、
浮上して水面で弾けるたびに、下から細いジェットが立って 10–100 μm の粒を
数 m/s で打ち上げる（炭酸飲料の上に漂うミストと同じ仕組み）。

参考:
- 巻き込み開始速度: Ervine et al. (1980), Zhu, Oguz & Prosperetti (2000, JFM 404),
  Kiger & Duncan (2012, Annu. Rev. Fluid Mech.)
- 巻き込み量: Bin (1993, Chem. Eng. Sci.) の相関を開始速度でゲート
- 気泡径: Hinze スケール d = 0.725 (σ/ρ)^0.6 ε^-0.4
- ジェット液滴の速度: Deike et al. (2018, PRF 3) V_d = 19 V_c (sqrt(La/La_c)-1)^-1/2 (1+2.2Bo)^-3/4
- ジェット液滴の径:   r_d/R = 0.6 (sqrt(La/La_c)-1)^(5/4) / La^(3/8)  (Gañán-Calvo 系の表式)
- 液滴の個数: Blanchard (1989) 海水で最大 6–7 個, R ≳ 1.5 mm で出なくなる
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .props import Fluid, G
from . import droplets


# ---------------------------------------------------------------- 巻き込み
V_E_SMOOTH = 3.0   # 乱れのない層流水柱の巻き込み開始速度 [m/s]（Zhu+2000: 乱れがなければ高速でも巻き込まない）
V_E_ROUGH = 0.8    # 乱れた（乱流）水柱の開始速度 [m/s]（Ervine 1980 / Kiger & Duncan 2012 の 0.8–1.1 m/s）
DELTA_REF = 0.08   # 表面乱れ（半径比）がこの程度で乱流相当とみなす


def onset_velocity(roughness: float, fluid: Fluid) -> float:
    """空気巻き込み開始速度 V_e [m/s]。表面の乱れが大きいほど低速で巻き込む。"""
    cap = math.sqrt(fluid.sigma / 0.0728)  # 表面張力が小さい（温水）ほど巻き込みやすい
    return cap * (V_E_ROUGH + (V_E_SMOOTH - V_E_ROUGH) * math.exp(-roughness / DELTA_REF))


def air_entrainment_ratio(U: float, d: float, H: float, d0: float, v0: float,
                          roughness: float, fluid: Fluid) -> float:
    """巻き込まれる空気の体積流量 / 水の体積流量。

    Bin (1993): Q_a/Q_w = 0.04 Fr0^0.28 (L/d0)^0.4 （乱流ジェット）
    を開始速度 V_e でゲートする: × S(U/V_e)，S は 1 付近で立ち上がる滑らかな関数。
    """
    Ve = onset_velocity(roughness, fluid)
    Fr0 = v0 ** 2 / (G * d0)
    base = 0.04 * Fr0 ** 0.28 * max(H / d0, 1.0) ** 0.4
    x = U / Ve
    # x<1 でも稀な大きな乱れで散発的に巻き込む（ロジスティックの裾）
    gate = (max(1.0 - 1.0 / x ** 2, 0.0)) ** 1.5 + 0.02 / (1.0 + math.exp(-(x - 1.0) / 0.06))
    return base * min(gate, 1.0)


def bubble_size_median(U: float, d: float, fluid: Fluid, C_eps: float = 0.05) -> float:
    """巻き込まれた気泡の中央半径 [m]（Hinze スケール）。"""
    eps = C_eps * U ** 3 / d
    dH = 0.725 * (fluid.sigma / fluid.rho) ** 0.6 * eps ** -0.4
    return float(np.clip(0.5 * dH, 30e-6, 2e-3))


def bubble_distribution(R_med: float, n: int = 64, sigma_g: float = 0.6):
    """対数正規の気泡半径分布（体積重み）。返り値: (半径, 体積分率)。"""
    lnR = np.linspace(math.log(R_med) - 3 * sigma_g, math.log(R_med) + 3 * sigma_g, n)
    w = np.exp(-0.5 * ((lnR - math.log(R_med)) / sigma_g) ** 2)
    w /= w.sum()
    return np.exp(lnR), w


# ---------------------------------------------------------------- 破裂ジェット液滴
LA_C = 500.0  # これより小さい Laplace 数ではジェット液滴が出ない（Deike+2018, Berny+2020）


def jet_drops(R, fluid: Fluid):
    """半径 R の気泡が弾けたときに出るジェット液滴（先頭液滴）の (半径, 速度, 個数)。"""
    R = np.asarray(R, float)
    La = R / fluid.visco_capillary_length
    Bo = fluid.bond(R)
    Vc = np.sqrt(fluid.sigma / (fluid.rho * R))
    s = np.sqrt(np.maximum(La / LA_C, 1.0 + 1e-9)) - 1.0
    Vd = 19.0 * Vc * s ** -0.5 * (1.0 + 2.2 * Bo) ** -0.75
    rd = R * 0.6 * s ** 1.25 / La ** 0.375
    # 液滴の数: 小さい気泡で最大 7 個, Bo が大きく（R ≳ 1.5 mm）なると 0
    n = np.clip(np.round(7.0 * (1.0 - Bo / 0.30)), 0, 7)
    n = np.where(La < LA_C, 0, n)
    return rd, Vd, n


@dataclass
class SprayResult:
    air_ratio: float          # 巻き込み空気 / 水 (体積比)
    bubble_rate: float        # 気泡の生成率 [個/s]
    R_med: float              # 気泡の中央半径 [m]
    drop_rate: float          # ジェット液滴の発生率 [個/s]
    drop_mass_rate: float     # ジェット液滴の質量発生率 [kg/s]
    h_drop_median: float      # 飛沫の到達高さ（中央値）[m]
    h_drop_p99: float         # 飛沫の到達高さ（上位 1%）[m]
    heights: np.ndarray       # 到達高さのサンプル [m]
    weights: np.ndarray       # その頻度 [個/s]
    radii: np.ndarray         # その粒半径 [m]

    def rate_above(self, h: float) -> float:
        """到達高さが h を超える飛沫の発生率 [個/s]。"""
        return float(self.weights[self.heights > h].sum())

    def mass_rate_above(self, h: float, rho: float = 998.0) -> float:
        m = rho * 4.0 / 3.0 * math.pi * self.radii ** 3
        return float((self.weights * m)[self.heights > h].sum())


def spray_from_bubbles(Q_air: float, R_med: float, fluid: Fluid) -> SprayResult:
    """空気の巻き込み率 Q_air [m^3/s] から、気泡破裂による飛沫の統計を出す。"""
    R, wvol = bubble_distribution(R_med)
    Vb = 4.0 / 3.0 * math.pi * R ** 3
    nb = Q_air * wvol / Vb                    # 各ビンの気泡生成率 [個/s]
    rd, Vd, ndrops = jet_drops(R, fluid)
    hs, ws, rs = [], [], []
    for k in range(7):
        # k 番目の液滴: 径は ~8% ずつ大きく、速度は ~0.6 倍ずつ遅い（Blanchard 1989, Berny 2021 の傾向）
        rk = rd * 1.08 ** k
        vk = Vd * 0.6 ** k
        has = ndrops > k
        hk = droplets.max_height(rk, vk, rho_l=fluid.rho)
        hs.append(np.where(has, hk, 0.0))
        ws.append(np.where(has, nb, 0.0))
        rs.append(rk)
    hs = np.concatenate(hs); ws = np.concatenate(ws); rs = np.concatenate(rs)
    drop_rate = float(ws.sum())
    mass_rate = float((ws * fluid.rho * 4 / 3 * math.pi * rs ** 3).sum())
    if drop_rate > 0:
        order = np.argsort(hs)
        cdf = np.cumsum(ws[order]) / drop_rate
        h50 = float(hs[order][np.searchsorted(cdf, 0.5)])
        h99 = float(hs[order][min(np.searchsorted(cdf, 0.99), len(cdf) - 1)])
    else:
        h50 = h99 = 0.0
    return SprayResult(air_ratio=float("nan"), bubble_rate=float(nb.sum()), R_med=R_med,
                       drop_rate=drop_rate, drop_mass_rate=mass_rate,
                       h_drop_median=h50, h_drop_p99=h99, heights=hs, weights=ws, radii=rs)
