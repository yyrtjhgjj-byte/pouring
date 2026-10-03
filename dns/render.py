"""DNS スナップショットの可視化（軸で折り返した断面図）と動画化。

使い方:
  python3 dns/render.py dns/runs/startup_H10cm            # 全フレーム → frames/*.png → movie.mp4
  python3 dns/render.py dns/runs/startup_H10cm --times 0.01 0.03 0.06 --montage out.png
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["IPAGothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
import snap  # noqa: E402

WATER = LinearSegmentedColormap.from_list("water", ["#0b1320", "#123a63", "#1f6fb2", "#7cc4ff"])
SPEED = LinearSegmentedColormap.from_list("speed", ["#0d2847", "#1565c0", "#4fc3f7", "#fff59d", "#ff7043"])


def frame(s: snap.Snap, ax, meta: dict | None = None, vmax: float = 2.5, zoom=None, title=True):
    """左半分: 速度の大きさ（水の中だけ）, 右半分: 体積率。中央が対称軸。"""
    c = s.mirrored("c")
    spd = np.hypot(s.mirrored("u"), s.mirrored("w"))
    nz, n2 = c.shape
    ext = [-s.nr * s.dx * 1e3, s.nr * s.dx * 1e3, 0, nz * s.dx * 1e3]
    half = n2 // 2
    img = np.zeros((nz, n2, 4))
    # 右: 水（体積率）
    rgba_c = WATER(np.clip(c, 0, 1))
    # 左: 速度（水の中）
    rgba_v = SPEED(np.clip(spd / vmax, 0, 1))
    air = np.array([0.96, 0.97, 0.99, 1.0])
    left = rgba_v * c[..., None] + air * (1 - c[..., None])
    right = rgba_c * c[..., None] + air * (1 - c[..., None])
    img[:, :half] = left[:, :half]
    img[:, half:] = right[:, half:]
    ax.imshow(np.clip(img, 0, 1), extent=ext, interpolation="bilinear")
    r = np.concatenate([-s.r[::-1], s.r]) * 1e3
    z = s.z * 1e3
    ax.contour(r, z, c[::-1, :], levels=[0.5], colors="#0a2540", linewidths=0.5)
    if meta and "pool" in meta:
        ax.axhline(meta["pool"] * 1e3, color="#999", lw=0.4, ls=":")
    if zoom:
        ax.set_xlim(zoom[0], zoom[1]); ax.set_ylim(zoom[2], zoom[3])
    ax.set_xlabel("r [mm]"); ax.set_ylabel("z [mm]")
    if title:
        ax.set_title(f"t = {s.t*1e3:6.1f} ms")
    ax.set_aspect("equal")


def load_meta(d: Path) -> dict:
    m = {}
    if (d / "meta.json").exists():
        m = json.load(open(d / "meta.json"))
    prm = {}
    if (d / "params.txt").exists():
        for line in open(d / "params.txt"):
            k, v = line.split()[:2]
            prm[k] = v
    if "pool_depth" in prm:
        m["pool"] = float(prm["pool_depth"])
    return m


def render_movie(d: Path, every: int = 1, fps: int = 30, vmax: float = 2.5, zoom=None, dpi=110):
    files = snap.list_snaps(str(d))[::every]
    meta = load_meta(d)
    out = d / "frames"
    out.mkdir(exist_ok=True)
    for k, f in enumerate(files):
        s = snap.read(f)
        fig, ax = plt.subplots(figsize=(6.4, 7.2))
        frame(s, ax, meta, vmax=vmax, zoom=zoom)
        sm = plt.cm.ScalarMappable(cmap=SPEED, norm=plt.Normalize(0, vmax))
        cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("水の速さ [m/s]（左半分）", fontsize=8)
        fig.tight_layout()
        fig.savefig(out / f"f_{k:05d}.png", dpi=dpi)
        plt.close(fig)
    mp4 = d / "movie.mp4"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(out / "f_%05d.png"),
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", "-c:v", "libx264",
                    "-crf", "23", str(mp4)], check=True)
    return mp4


def montage(d: Path, times, out: Path, zoom=None, vmax=2.5, ncol=None, suptitle=None):
    files = snap.list_snaps(str(d))
    ts = np.array([snap.read(f).t for f in files])
    meta = load_meta(d)
    n = len(times)
    ncol = ncol or n
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 3.4 * nrow), squeeze=False)
    for k, t in enumerate(times):
        s = snap.read(files[int(np.argmin(np.abs(ts - t)))])
        ax = axs[k // ncol][k % ncol]
        frame(s, ax, meta, vmax=vmax, zoom=zoom)
        if k % ncol:
            ax.set_ylabel("")
    for k in range(n, nrow * ncol):
        axs[k // ncol][k % ncol].axis("off")
    if suptitle:
        fig.suptitle(suptitle)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--times", type=float, nargs="*")
    ap.add_argument("--montage", default=None)
    ap.add_argument("--zoom", type=float, nargs=4, default=None)
    ap.add_argument("--vmax", type=float, default=2.5)
    a = ap.parse_args()
    d = Path(a.run)
    if a.times:
        print(montage(d, a.times, Path(a.montage or d / "montage.png"), zoom=a.zoom, vmax=a.vmax))
    else:
        print(render_movie(d, every=a.every, zoom=a.zoom, vmax=a.vmax))
