"""axi2phase のスナップショット読み込み・解析ユーティリティ。"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd  # noqa: F401  (diag の読み込みで使う)
from scipy import ndimage


@dataclass
class Snap:
    t: float
    dx: float
    c: np.ndarray   # [nr, nz] 体積率
    u: np.ndarray   # セル中心の r 速度
    w: np.ndarray   # セル中心の z 速度
    p: np.ndarray

    @property
    def nr(self):
        return self.c.shape[0]

    @property
    def nz(self):
        return self.c.shape[1]

    @property
    def r(self):
        return (np.arange(self.nr) + 0.5) * self.dx

    @property
    def z(self):
        return (np.arange(self.nz) + 0.5) * self.dx

    def mirrored(self, field: str = "c") -> np.ndarray:
        """軸で折り返した [nz, 2nr] 画像（上が +z になるよう反転済み）。"""
        f = getattr(self, field)
        full = np.concatenate([f[::-1, :], f], axis=0)  # [2nr, nz]
        return full.T[::-1, :]


def read(fn: str) -> Snap:
    with open(fn, "rb") as f:
        nr, nz = np.frombuffer(f.read(8), dtype=np.int32)
        dx, t = np.frombuffer(f.read(16), dtype=np.float64)
        data = np.frombuffer(f.read(), dtype=np.float32).reshape(4, nr, nz)
    return Snap(t=float(t), dx=float(dx), c=data[0], u=data[1], w=data[2], p=data[3])


def list_snaps(d: str) -> list[str]:
    return sorted(glob.glob(os.path.join(d, "snap_*.bin")))


def read_diag(d: str):
    return pd.read_csv(os.path.join(d, "diag.csv"))


def liquid_structures(s: Snap, z_min: float = 0.0, thresh: float = 0.5):
    """液体の連結成分（滴・柱）を抽出し、z_min より上にある部分の統計を返す。

    返り値: list of dict(volume, z_top, z_bottom, z_cm, r_cm, w_cm, w_top, eq_radius)
    体積は軸対称 2πr 重み。
    """
    mask = s.c > thresh
    lab, n = ndimage.label(mask)
    out = []
    if n == 0:
        return out
    R = s.r[:, None] * np.ones((1, s.nz))
    Z = np.ones((s.nr, 1)) * s.z[None, :]
    vol_w = 2 * np.pi * R * s.c * s.dx ** 2
    idx = np.arange(1, n + 1)
    vol = ndimage.sum(vol_w, lab, idx)
    zc = ndimage.sum(vol_w * Z, lab, idx) / np.maximum(vol, 1e-30)
    rcm = ndimage.sum(vol_w * R, lab, idx) / np.maximum(vol, 1e-30)
    wc = ndimage.sum(vol_w * s.w, lab, idx) / np.maximum(vol, 1e-30)
    objs = ndimage.find_objects(lab)
    for k, sl in enumerate(objs):
        zt = (sl[1].stop) * s.dx
        zb = (sl[1].start) * s.dx
        if zt < z_min:
            continue
        sub = lab[sl] == (k + 1)
        # 最上部の速度
        jtop = sl[1].stop - 1
        col = (lab[:, jtop] == (k + 1))
        wtop = float(np.mean(s.w[col, jtop])) if col.any() else 0.0
        out.append({"label": k + 1, "volume": float(vol[k]), "z_top": zt, "z_bottom": zb,
                    "z_cm": float(zc[k]), "r_cm": float(rcm[k]), "w_cm": float(wc[k]),
                    "w_top": wtop, "r_max": sl[0].stop * s.dx,
                    "eq_radius": float((3 * vol[k] / (4 * np.pi)) ** (1 / 3))})
    return out
