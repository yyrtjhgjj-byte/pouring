"""軸対称 DNS の結果を 3D に回転させて、手前半分を切り取った「断面つき 3D 画像」を描く。

- 界面（c = 0.5）を軸まわりに回転させた面を、2D の符号付き距離場で sphere tracing
- 切断面（y = 0）は体積率をそのまま塗る → 水中の空洞や気泡が見える
- 水面は Fresnel 反射つきの簡易シェーディング、奥半分のコップ（ガラス）も薄く描く

  python3 dns/render3d.py dns/runs/startup_H10cm --t 0.06 --out frame.png
  python3 dns/render3d.py dns/runs/startup_H10cm --movie  → movie3d.mp4
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from numba import njit, prange

sys.path.insert(0, str(Path(__file__).resolve().parent))
import snap  # noqa: E402


def sdf2d(c: np.ndarray, dx: float) -> np.ndarray:
    """子午面 (r,z) の符号付き距離（液体内が負）[m]。"""
    liq = c > 0.5
    d_out = ndimage.distance_transform_edt(~liq) * dx
    d_in = ndimage.distance_transform_edt(liq) * dx
    # セル内の界面位置で 0.5 セル補正
    return np.where(liq, -(d_in - 0.5 * dx), d_out - 0.5 * dx)


@njit(cache=True, fastmath=True)
def bilin(f, dx, r, z):
    nr, nz = f.shape
    x = r / dx - 0.5
    y = z / dx - 0.5
    if x < 0.0:
        x = 0.0
    if y < 0.0:
        y = 0.0
    if x > nr - 1.001:
        x = nr - 1.001
    if y > nz - 1.001:
        y = nz - 1.001
    i = int(x)
    j = int(y)
    tx = x - i
    ty = y - j
    return (f[i, j] * (1 - tx) * (1 - ty) + f[i + 1, j] * tx * (1 - ty)
            + f[i, j + 1] * (1 - tx) * ty + f[i + 1, j + 1] * tx * ty)


@njit(cache=True, fastmath=True)
def env(dx_, dy_, dz_):
    # スタジオ風の環境: 上が明るく、下は暗い。左上に大きなソフトボックス
    t = 0.5 * (dz_ + 1.0)
    r = 0.06 + 0.30 * t
    g = 0.09 + 0.36 * t
    b = 0.14 + 0.45 * t
    s = max(0.0, -0.55 * dx_ + 0.25 * dy_ + 0.8 * dz_)
    s = s ** 24 * 2.2
    return r + s, g + s, b + s


@njit(parallel=True, cache=True, fastmath=True)
def render(c, cs, sd, dx, Lr, Lz, z_pool, W, H, cam, look, fov, out):
    # カメラ基底
    fwd = look - cam
    fwd /= np.sqrt((fwd ** 2).sum())
    up0 = np.array([0.0, 0.0, 1.0])
    right = np.cross(fwd, up0)
    right /= np.sqrt((right ** 2).sum())
    up = np.cross(right, fwd)
    tanf = np.tan(fov * 0.5)
    R_cup = Lr
    for py in prange(H):
        for px in range(W):
            u = (2.0 * (px + 0.5) / W - 1.0) * tanf * W / H
            v = (1.0 - 2.0 * (py + 0.5) / H) * tanf
            d = fwd + u * right + v * up
            d /= np.sqrt((d ** 2).sum())
            # 背景
            br, bg, bb = 0.05 + 0.03 * v, 0.09 + 0.05 * v, 0.14 + 0.07 * v
            col_r, col_g, col_b = br, bg, bb
            # ドメイン（円柱 r<Lr, 0<z<Lz）と y>=0 の半分の交差区間
            ox, oy, oz = cam[0], cam[1], cam[2]
            dx_, dy_, dz_ = d[0], d[1], d[2]
            # 円柱との交点
            a = dx_ * dx_ + dy_ * dy_
            bq = 2 * (ox * dx_ + oy * dy_)
            cq = ox * ox + oy * oy - R_cup * R_cup
            disc = bq * bq - 4 * a * cq
            if disc < 0 or a < 1e-12:
                out[py, px, 0] = col_r; out[py, px, 1] = col_g; out[py, px, 2] = col_b
                continue
            sq = np.sqrt(disc)
            t0 = (-bq - sq) / (2 * a)
            t1 = (-bq + sq) / (2 * a)
            # z 範囲
            if abs(dz_) > 1e-12:
                tz0 = (0.0 - oz) / dz_
                tz1 = (Lz - oz) / dz_
                if tz0 > tz1:
                    tz0, tz1 = tz1, tz0
                t0 = max(t0, tz0)
                t1 = min(t1, tz1)
            # y >= 0 の半空間（手前 y<0 を切り取る）
            if abs(dy_) > 1e-12:
                ty = -oy / dy_
                if dy_ > 0:
                    t0 = max(t0, ty)
                else:
                    t1 = min(t1, ty)
            elif oy < 0:
                t1 = -1.0
            hit = False
            if t1 > t0 and t1 > 0:
                t = max(t0, 0.0) + 1e-7
                # 入口が切断面 (y=0) 上で、そこが水なら断面を塗る
                pxx = ox + t * dx_; pzz = oz + t * dz_
                pyy = oy + t * dy_
                rr = np.sqrt(pxx * pxx + pyy * pyy)
                on_cut = abs(pyy) < 1e-6
                if on_cut and bilin(c, dx, rr, pzz) > 0.5:
                    # 断面: 体積率と深さで濃淡
                    cf = bilin(c, dx, rr, pzz)
                    depth = max(0.0, (z_pool - pzz) / z_pool)
                    col_r = 0.10 + 0.10 * (1 - depth)
                    col_g = 0.40 + 0.18 * (1 - depth)
                    col_b = 0.72 + 0.12 * (1 - depth)
                    # 断面の縁（界面近く）を明るく
                    sdv = abs(bilin(sd, dx, rr, pzz))
                    edge = max(0.0, 1.0 - sdv / (2.5 * dx))
                    col_r += 0.45 * edge; col_g += 0.45 * edge; col_b += 0.28 * edge
                    hit = True
                else:
                    for it in range(600):
                        pxx = ox + t * dx_; pyy = oy + t * dy_; pzz = oz + t * dz_
                        rr = np.sqrt(pxx * pxx + pyy * pyy)
                        s = bilin(sd, dx, rr, pzz)
                        if s < 0.25 * dx:
                            # 二分法で界面に寄せる
                            ta = t - max(abs(s), dx); tb = t
                            for k in range(10):
                                tm = 0.5 * (ta + tb)
                                qx = ox + tm * dx_; qy = oy + tm * dy_; qz = oz + tm * dz_
                                if bilin(c, dx, np.sqrt(qx * qx + qy * qy), qz) > 0.5:
                                    tb = tm
                                else:
                                    ta = tm
                            t = tb
                            pxx = ox + t * dx_; pyy = oy + t * dy_; pzz = oz + t * dz_
                            rr = np.sqrt(pxx * pxx + pyy * pyy) + 1e-12
                            # 法線: 子午面の勾配を 3D に
                            e = 1.0 * dx
                            gr = bilin(cs, dx, rr + e, pzz) - bilin(cs, dx, max(rr - e, 0.0), pzz)
                            gz = bilin(cs, dx, rr, pzz + e) - bilin(cs, dx, rr, pzz - e)
                            nx = -gr * pxx / rr; ny = -gr * pyy / rr; nz = -gz
                            nn = np.sqrt(nx * nx + ny * ny + nz * nz) + 1e-12
                            nx /= nn; ny /= nn; nz /= nn
                            ndv = -(nx * dx_ + ny * dy_ + nz * dz_)
                            if ndv < 0:
                                nx = -nx; ny = -ny; nz = -nz; ndv = -ndv
                            # 反射ベクトル
                            rx = dx_ + 2 * ndv * nx; ry = dy_ + 2 * ndv * ny; rz = dz_ + 2 * ndv * nz
                            er, eg, eb = env(rx, ry, rz)
                            F = 0.02 + 0.98 * (1 - ndv) ** 5
                            lam = max(0.0, -0.45 * nx + 0.2 * ny + 0.87 * nz)
                            base_r = 0.06 + 0.22 * lam
                            base_g = 0.30 + 0.35 * lam
                            base_b = 0.55 + 0.35 * lam
                            col_r = base_r * (1 - F) + er * F
                            col_g = base_g * (1 - F) + eg * F
                            col_b = base_b * (1 - F) + eb * F
                            hit = True
                            break
                        t += max(s, 0.3 * dx)
                        if t > t1:
                            break
            # 奥半分のガラスのコップ（r = Lr, y > 0）
            if t1 > 0 and disc >= 0:
                tg = (-bq + sq) / (2 * a)
                gz_ = oz + tg * dz_
                gy_ = oy + tg * dy_
                if gy_ > 0 and gz_ > 0 and gz_ < Lz * 0.97:
                    gl = 0.10
                    col_r = col_r * (1 - gl) + 0.75 * gl
                    col_g = col_g * (1 - gl) + 0.85 * gl
                    col_b = col_b * (1 - gl) + 0.95 * gl
                    # ガラスの縁の光
                    gx_ = ox + tg * dx_
                    fr = abs(gx_) / R_cup
                    if fr > 0.97:
                        col_r += 0.25; col_g += 0.28; col_b += 0.3
            out[py, px, 0] = col_r
            out[py, px, 1] = col_g
            out[py, px, 2] = col_b


def render_snapshot(s: snap.Snap, z_pool: float, W=960, H=720, elev=20.0, azim=-90.0,
                    dist=0.16, fov=34.0, zoom_center=None):
    c = np.ascontiguousarray(s.c.astype(np.float64))
    sd = np.ascontiguousarray(sdf2d(s.c, s.dx).astype(np.float64))
    Lr, Lz = s.nr * s.dx, s.nz * s.dx
    zc = zoom_center if zoom_center is not None else 0.55 * Lz
    look = np.array([0.0, 0.0, zc])
    el, az = np.radians(elev), np.radians(azim)
    cam = look + dist * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    cs = np.ascontiguousarray(ndimage.gaussian_filter(s.c.astype(np.float64), 1.2, mode="nearest"))
    out = np.zeros((H, W, 3))
    render(c, cs, sd, s.dx, Lr, Lz, z_pool, W, H, cam, look, np.radians(fov), out)
    img = np.clip(out, 0, 1) ** (1 / 1.2)
    return (img * 255).astype(np.uint8)


def save_png(img, path, label=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    h, w, _ = img.shape
    fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
    if label:
        ax.text(0.02, 0.96, label, transform=ax.transAxes, color="#cfe3f5", fontsize=13,
                family="IPAGothic", va="top")
    fig.savefig(path, dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    import json
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--t", type=float, nargs="*", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--movie", action="store_true")
    ap.add_argument("--W", type=int, default=960)
    ap.add_argument("--H", type=int, default=720)
    ap.add_argument("--elev", type=float, default=18.0)
    ap.add_argument("--dist", type=float, default=0.15)
    ap.add_argument("--zc", type=float, default=None)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    d = Path(a.run)
    prm = {l.split()[0]: l.split()[1] for l in open(d / "params.txt")}
    z_pool = float(prm.get("pool_depth", 0.04))
    files = snap.list_snaps(str(d))
    ts = np.array([snap.read(f).t for f in files]) if not a.movie else None
    if a.movie:
        out = d / "frames3d"; out.mkdir(exist_ok=True)
        for k, f in enumerate(files):
            s = snap.read(f)
            img = render_snapshot(s, z_pool, a.W, a.H, elev=a.elev, dist=a.dist, zoom_center=a.zc)
            save_png(img, out / f"f_{k:05d}.png", f"{a.label}  t = {s.t*1e3:5.1f} ms")
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "30", "-i", str(out / "f_%05d.png"),
                        "-pix_fmt", "yuv420p", "-c:v", "libx264", "-crf", "22", str(d / "movie3d.mp4")], check=True)
        print(d / "movie3d.mp4")
    else:
        for t in a.t:
            k = int(np.argmin(np.abs(ts - t)))
            s = snap.read(files[k])
            img = render_snapshot(s, z_pool, a.W, a.H, elev=a.elev, dist=a.dist, zoom_center=a.zc)
            p = a.out or str(d / f"r3d_{t*1000:.0f}ms.png")
            if len(a.t) > 1:
                p = p.replace(".png", f"_{t*1000:.0f}ms.png")
            save_png(img, p, f"{a.label}  t = {s.t*1e3:5.1f} ms")
            print(p)
