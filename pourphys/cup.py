"""コップの幾何学：傾けたときの水面・こぼれ限界・水流の着地点・ノズルとの隙間。

座標系（世界系）: ノズル先端が原点、z は上向き、水流は x=0, y=0 を真下に落ちる。
コップ系 (X', Y', Z'): Z' がコップの軸（底 Z'=0, 縁 Z'=h）。
傾き θ > 0 で、コップの口が −x 側（左）に倒れる。このとき
  - 左側の内壁が「床」っぽく上を向き、水流を受け止められる
  - 左側の縁が一番低くなり、こぼれるのはここから
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class Cup:
    r_bottom: float = 0.030   # 底の内半径 [m]
    r_top: float = 0.037      # 縁の内半径 [m]
    height: float = 0.095     # 内側の高さ [m]（この既定値で満水 ~335 mL のタンブラー）

    @property
    def taper(self) -> float:
        """dr/dZ'（縁に向かって広がる割合）。"""
        return (self.r_top - self.r_bottom) / self.height

    def radius(self, Zp):
        return self.r_bottom + self.taper * np.asarray(Zp)

    @property
    def volume(self) -> float:
        h, a, b = self.height, self.r_bottom, self.r_top
        return math.pi * h / 3.0 * (a * a + a * b + b * b)

    def upright_level(self, V: float) -> float:
        """直立時に体積 V を入れたときの水深 [m]。"""
        lo, hi = 0.0, self.height
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            a, b = self.r_bottom, self.radius(mid)
            if math.pi * mid / 3.0 * (a * a + a * b + b * b) < V:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def _segment_area(r, d):
    """半径 r の円のうち、弦 X' < d の部分の面積（配列可）。"""
    r = np.asarray(r, float)
    d = np.clip(np.asarray(d, float), -r, r)
    return r * r * np.arccos(-d / r) + d * np.sqrt(np.maximum(r * r - d * d, 0.0))


@dataclass
class Pose:
    """コップの姿勢: 傾き θ [rad] と、底面中心の世界座標 (xc, zc)。"""

    cup: Cup
    theta: float
    xc: float = 0.0
    zc: float = -0.15

    # --- 座標変換 ---
    def to_world(self, Xp, Zp, Yp=0.0):
        c, s = math.cos(self.theta), math.sin(self.theta)
        Xp = np.asarray(Xp, float); Zp = np.asarray(Zp, float)
        return Xp * c - Zp * s + self.xc, np.asarray(Yp, float) + 0 * Xp, Xp * s + Zp * c + self.zc

    # --- 水 ---
    def water_volume_below(self, zw: float, n: int = 400) -> float:
        """世界系で水平面 z = zw より下にあるコップ内の体積 [m^3]。"""
        cup = self.cup
        Zp = (np.arange(n) + 0.5) * cup.height / n
        r = cup.radius(Zp)
        c, s = math.cos(self.theta), math.sin(self.theta)
        zc_slice = self.zc + Zp * c
        if abs(s) < 1e-12:
            A = np.where(zc_slice < zw, math.pi * r * r, 0.0)
        else:
            X_lim = (zw - zc_slice) / s
            A = _segment_area(r, X_lim) if s > 0 else math.pi * r * r - _segment_area(r, X_lim)
        return float(A.sum() * cup.height / n)

    def water_level(self, V: float) -> float:
        """体積 V の水を入れたときの水面の世界 z。"""
        _, _, zb = self.to_world(np.array([-self.cup.r_bottom, self.cup.r_bottom]), np.zeros(2))
        lo = float(np.min(zb)) - 1e-3
        hi = self.zc + self.cup.height + self.cup.r_top + 0.05
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self.water_volume_below(mid) < V:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def rim_points(self, n: int = 180):
        psi = np.linspace(0, 2 * math.pi, n, endpoint=False)
        Xp = self.cup.r_top * np.cos(psi)
        Yp = self.cup.r_top * np.sin(psi)
        return self.to_world(Xp, np.full(n, self.cup.height), Yp)

    def lowest_rim_z(self) -> float:
        return float(np.min(self.rim_points()[2]))

    def highest_rim_z(self) -> float:
        return float(np.max(self.rim_points()[2]))

    def max_volume(self) -> float:
        """この傾きでこぼれずに入る最大体積。"""
        return self.water_volume_below(self.lowest_rim_z())

    # --- 左（床側）の内壁 ---
    @property
    def wall_angle(self) -> float:
        """床側内壁の、鉛直からの傾き φ_w [rad] = θ + テーパー角。"""
        return self.theta + math.atan(self.cup.taper)

    def left_wall_point(self, Zp):
        return self.to_world(-self.cup.radius(Zp), Zp)

    def left_wall_Zp_at_x(self, x: float) -> float:
        """床側内壁が世界 x を横切るコップ系 Z'。"""
        c, s, k = math.cos(self.theta), math.sin(self.theta), self.cup.taper
        return (self.xc - x - self.cup.r_bottom * c) / (k * c + s)

    def left_wall_Zp_at_z(self, z: float) -> float:
        c, s, k = math.cos(self.theta), math.sin(self.theta), self.cup.taper
        # z = -r(Z') s + Z' c + zc
        return (z - self.zc + self.cup.r_bottom * s) / (c - k * s)

    # --- ノズルとの隙間 ---
    def nozzle_clearance(self, r_nozzle: float = 0.007, nozzle_len: float = 0.05) -> float:
        """縁とノズル（半径 r_nozzle の縦の円柱, 先端 z=0）の最短距離 [m]。負なら干渉。"""
        x, y, z = self.rim_points(360)
        rho = np.hypot(x, y)
        out_r = np.maximum(rho - r_nozzle, 0.0)
        below = np.maximum(-z, 0.0)
        above = np.maximum(z - nozzle_len, 0.0)
        d = np.sqrt(out_r ** 2 + below ** 2 + above ** 2)
        inside = (rho < r_nozzle) & (z >= 0) & (z <= nozzle_len)
        d = np.where(inside, -np.minimum(r_nozzle - rho, z), d)
        # 縁がノズル先端より上にある＝ノズルがコップの口の中に入っている → 壁までの水平距離で評価
        enclosing = (z > 0).any() and (rho > r_nozzle).all()
        if enclosing:
            return float(np.min(rho - r_nozzle))
        return float(np.min(d))


def place_cup(cup: Cup, theta: float, V: float, gap: float, aim: str = "center",
              wall_offset: float = 0.01, x_aim: float | None = None) -> tuple[Pose, dict]:
    """注ぎ方（傾き・隙間・狙い）からコップの姿勢を決める。

    gap: ノズル先端から、縁の一番高い点までの鉛直距離 [m]（正 = 縁がノズルより下）
    aim: "center"（水面の中央に落とす）, "wall"（床側の内壁の、水面から wall_offset 上に当てる）
    """
    pose = Pose(cup, theta, 0.0, 0.0)
    zw_local = pose.water_level(V)  # 姿勢 (xc=0,zc=0) での水面
    info = {"spill": V > pose.max_volume() + 1e-12}
    if aim == "wall":
        Zp_w = pose.left_wall_Zp_at_z(zw_local)
        dZ = wall_offset / math.cos(pose.wall_angle) / math.sqrt(1 + cup.taper ** 2)
        Zp_hit = min(Zp_w + dZ, cup.height * 0.98)
        xw, _, _ = pose.left_wall_point(Zp_hit)
        xc = -float(xw)  # 水流 (x=0) がそこを通るよう平行移動
    else:
        # 水面の（x 方向の）中点を狙う: 床側内壁と反対側の壁の水面位置の中点
        Zp_w = pose.left_wall_Zp_at_z(zw_local)
        xl, _, _ = pose.left_wall_point(Zp_w)
        c, s, k = math.cos(pose.theta), math.sin(pose.theta), cup.taper
        # 右壁: X' = +r(Z'),  z = r s + Z' c → Z' = (z - r_b s)/(c + k s)
        Zp_r = (zw_local - cup.r_bottom * s) / (c + k * s)
        Zp_r = min(max(Zp_r, 0.0), cup.height)
        xr, _, _ = pose.to_world(cup.radius(Zp_r), Zp_r)
        xc = -0.5 * (float(xl) + float(xr))
    if x_aim is not None:
        xc += x_aim
    pose = Pose(cup, theta, xc, 0.0)
    zc = -gap - pose.highest_rim_z()
    pose = Pose(cup, theta, xc, zc)
    info["z_water"] = pose.water_level(V)
    info["lowest_rim"] = pose.lowest_rim_z()
    info["highest_rim"] = pose.highest_rim_z()
    return pose, info
