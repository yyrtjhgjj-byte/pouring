"""注ぎ始めの過渡現象：「ジェット頭部（先頭の水塊）」の形成。

レバーを開けた直後、吐出速度は 0 から v0 まで立ち上がる。
後から出た速い水が、先に出た遅い水に追いつくと合体して塊（頭部）になる。
これを「粘着粒子（sticky particle）」＝ 付着型 Burgers 方程式のラグランジュ解法で解く。

  各水塊 i: 放出時刻 t0_i, 初速 v_i, 質量 m_i
  重力は全員に等しく効くので、自由落下系では等速直線運動
  追突したら運動量保存で合体（完全非弾性）

さらに、表面張力で先端が丸まる「毛管バルブ」（直径 ~1.6 d）も下限として考慮する。

この頭部が水面に突っ込むと、クレーターが潰れて上向きの細いジェット
（Worthington ジェット＝「水滴の柱」）が立つ。Michon et al. (2017) によると
液滴衝突では Fr = U²/(gD) ≳ 60 でジェットが現れ、Fr ≳ 90 で衝突速度の 3.5–4 倍に達する。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .props import G
from .server import Server


@dataclass
class HeadResult:
    H: float            # 落下距離 [m]
    mass: float         # 頭部の質量 [kg]
    U: float            # 着水時の頭部速度 [m/s]
    D: float            # 頭部の等価直径 [m]
    t_arrive: float     # 開栓から着水までの時間 [s]
    kinematic: bool     # 追突合体で頭部ができたか（False なら毛管バルブのみ）
    n_merged: int       # 合体した水塊の数
    drip: bool          # 立ち上がりが遅く、最初は滴下（ポタポタ）になるか

    def froude(self) -> float:
        return self.U ** 2 / (G * self.D)


def simulate_head(server: Server, H: float, n_parcels: int = 1500,
                  t_release: float | None = None, dt: float = 2e-5) -> HeadResult:
    """粘着粒子モデルで、落下距離 H に頭部が届いた瞬間の状態を求める。"""
    fl = server.fluid
    tau = max(server.tau_valve, 1e-4)
    if t_release is None:
        # 着水までに放出される分だけ追えばよい
        t_release = min(tau * 3 + 0.05, 0.6)
    t0 = (np.arange(n_parcels) + 0.5) * (t_release / n_parcels)
    v = server.exit_velocity(t0)
    m = fl.rho * server.flow(t0) * (t_release / n_parcels)

    # 滴下判定: 吐出口での We < 4 の間はジェットにならず滴になる
    We_exit = fl.rho * v ** 2 * server.d0 / fl.sigma
    drip = bool(np.any((We_exit < 4.0) & (t0 < tau)) and tau > 0.15)

    # 有効な（流量がある）水塊のみ
    ok = m > 1e-12
    t0, v, m = t0[ok], v[ok], m[ok]
    # 状態: 位置 z（下向き）, 速度 w。放出前の水塊は z=0 に待機
    z = np.zeros_like(v)
    w = v.copy()
    alive = np.ones_like(v, dtype=bool)
    nmerge = np.ones_like(v, dtype=int)
    released = np.zeros_like(v, dtype=bool)

    t = 0.0
    lead = 0
    while t < 2.0:
        t += dt
        act = released & alive
        z[act] += w[act] * dt + 0.5 * G * dt * dt
        w[act] += G * dt
        # このステップ内で放出された水塊は、放出時刻からの経過分だけ進めて登場させる
        newly = (~released) & (t0 <= t)
        if np.any(newly):
            s = t - t0[newly]
            z[newly] = v[newly] * s + 0.5 * G * s * s
            w[newly] = v[newly] + G * s
            released |= newly
            act = released & alive
        # 追突判定: 放出順（t0 の小さい順＝先頭側）に並んでいるので、
        # 後続 j が前方 i を追い越したら合体
        idx = np.flatnonzero(act)
        if idx.size >= 2:
            zi = z[idx]
            # 前（index 小）より後ろ（index 大）の z が大きければ追い越し
            cross = np.flatnonzero(zi[1:] >= zi[:-1])
            for c in cross[::-1]:
                i, j = idx[c], idx[c + 1]
                if not (alive[i] and alive[j]):
                    continue
                M = m[i] + m[j]
                w[i] = (m[i] * w[i] + m[j] * w[j]) / M
                z[i] = (m[i] * z[i] + m[j] * z[j]) / M
                m[i] = M
                nmerge[i] += nmerge[j]
                alive[j] = False
        live = np.flatnonzero(alive & released)
        if live.size == 0:
            continue
        lead = int(live[0])
        if z[lead] >= H:
            break

    # 先頭から「毛管バルブ 1 個分（直径 1.6 d の球）」の質量までを頭部とみなし、
    # その運動量平均速度で衝突させる。追突合体の塊がそれより重ければ塊全体を頭部とする。
    U_front = float(w[lead])
    d_local = math.sqrt(4 * server.Q / (math.pi * max(U_front, 1e-6)))
    m_bulb = fl.rho * math.pi / 6.0 * (1.6 * d_local) ** 3
    live = np.flatnonzero(alive & released)
    order = live[np.argsort(-z[live])]
    cm = np.cumsum(m[order])
    k = int(np.searchsorted(cm, m_bulb)) + 1
    sel = order[:max(k, 1)]
    M = float(m[sel].sum())
    # 各水塊が H に届いたときの速度（そこまで弾道で落ちると近似）
    w_at_H = np.sqrt(w[sel] ** 2 + 2.0 * G * np.maximum(H - z[sel], 0.0))
    U = float((m[sel] * w_at_H).sum() / M)
    kinematic = nmerge[lead] > max(3, 0.02 * len(m)) and m[lead] > 0.5 * m_bulb
    D_eq = (6.0 * M / (fl.rho * math.pi)) ** (1.0 / 3.0)
    # 頭部は最大でも水柱径の 2.5 倍程度の丸い塊にしかならない（それ以上は細長い）
    D = float(np.clip(D_eq, d_local, 2.5 * d_local))
    return HeadResult(H=H, mass=M, U=U, D=D, t_arrive=t,
                      kinematic=bool(kinematic), n_merged=int(nmerge[lead]), drip=drip)


def worthington_regime(Fr: float) -> str:
    """Michon et al. (2017) の液滴衝突レジーム分類（深いプール）。"""
    if Fr < 60.0:
        return "none"        # ジェットなし（またはごく低い盛り上がり）
    if Fr < 90.0:
        return "slow"        # 低速のジェット
    return "fast"            # 高速・細いジェット（衝突速度の 3.5–4 倍）


def worthington_jet_velocity(U: float, Fr: float) -> float:
    """Worthington ジェット先端速度の区分的なモデル [m/s]。

    文献値（Fr<60 なし, 60–90 で立ち上がり, ≳90 で 3.5–4U）を
    滑らかにつないだ経験的補間。DNS（dns/）の結果で後から係数を較正する。
    """
    if Fr < 45.0:
        # 盛り上がり（低い「こぶ」）程度：重力波的な戻り ~0.3U
        return 0.3 * U * (Fr / 45.0)
    if Fr < 60.0:
        return U * (0.3 + 0.4 * (Fr - 45.0) / 15.0)
    if Fr < 90.0:
        return U * (0.7 + 2.8 * ((Fr - 60.0) / 30.0) ** 2)
    return 3.75 * U
