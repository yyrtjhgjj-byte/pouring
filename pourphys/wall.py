"""傾けたコップの内壁に水流を当てたときの「壁面膜流れ」モデル。

ビールを泡立てずに注ぐときの、グラスを傾けて内壁に沿わせるやり方。

1. 斜め衝突: 水流は壁に当たると薄い膜になって広がる。非粘性なら速さの大きさは保存されるが、
   壁に垂直な成分は失われて横・前後方向の広がりに変わる（衝突で粒は飛ばない:
   濡れ壁への液滴衝突のスプラッシュ指標 K = We^0.5 Re^0.25 は ~100 で、閾値 ~2000 より遥かに小さい）。
2. 膜の流下: 壁に沿った重力加速と、粘性摩擦（半放物線速度分布 τ = 3μu/h）。
3. 水面への進入: 膜は壁に沿って斜めに滑り込む。水面に垂直な成分 u cos φ_w だけが
   「突っ込む」速度になり、膜の表面は水柱よりずっと滑らか（水柱の凸凹は壁に付いた時点で均される）。
   さらに膜の裏側は壁なので、水柱のように裏側から空気の空洞ができない。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .props import Fluid, G


@dataclass
class FilmResult:
    u_impact: float      # 壁衝突時の水流速度 [m/s]
    u_entry: float       # 水面進入時の膜速度 [m/s]
    h_entry: float       # 膜厚 [m]
    width: float         # 膜幅 [m]
    U_normal: float      # 水面に垂直な進入速度成分 [m/s]
    roughness: float     # 膜表面の相対的な乱れ
    path: float          # 壁を流れた距離 [m]
    splash_K: float      # 壁衝突のスプラッシュ指標（参考）
    offset_x: float      # 着水点の水平ずれ（壁衝突点からの）[m]


def wall_film(U: float, a: float, Q: float, phi_w: float, s_path: float,
              jet_roughness: float, fluid: Fluid, n: int = 400) -> FilmResult:
    """壁衝突から水面進入までの膜。

    U: 壁に当たる速度, a: 水柱半径, Q: 流量, phi_w: 壁の鉛直からの傾き [rad],
    s_path: 衝突点から水面までの壁沿い距離 [m]
    """
    # 濡れ壁への衝突スプラッシュ指標（壁に垂直な成分で評価）
    Un = U * math.sin(phi_w)
    d = 2 * a
    K = math.sqrt(max(fluid.weber(Un, d), 0.0)) * max(fluid.reynolds(Un, d), 1e-9) ** 0.25
    # 膜幅: かすめる衝突ほど細長く、正面衝突ほど横に広がる
    w = d * (1.0 + 3.0 * math.sin(phi_w))
    # 衝突で失う分: 垂直成分の運動エネルギーの一部は横方向の広がり＋渦で散逸（係数 0.5 と仮定）
    u = math.sqrt(max(U * U - 0.5 * Un * Un, 1e-8))
    gs = G * math.cos(phi_w)
    ds = s_path / n
    s = 0.0
    for _ in range(n):
        h = Q / (u * w)
        du = (gs - 3.0 * fluid.nu * u / (h * h)) / u * ds
        u = max(u + du, 0.05)
        s += ds
    h = Q / (u * w)
    rough = 0.25 * jet_roughness  # 壁に付くと水柱の凸凹はほぼ均される（仮定）
    return FilmResult(u_impact=U, u_entry=u, h_entry=h, width=w,
                      U_normal=u * math.cos(phi_w), roughness=rough, path=s_path,
                      splash_K=K, offset_x=s_path * math.sin(phi_w))
