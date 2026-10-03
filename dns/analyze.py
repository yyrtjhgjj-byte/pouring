"""DNS の結果解析：水面より上に飛んだ水滴、立ち上がる水の柱、巻き込まれた空気。

  python3 dns/analyze.py dns/runs/startup_H10cm   → analysis.json と timeseries.csv

軸対称なので、軸から離れた位置の「滴」は実際にはドーナツ状のリング（3D では多数の滴に分かれる）。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import snap  # noqa: E402
from pourphys import droplets  # noqa: E402

G = 9.80665


def load_params(d: Path) -> dict:
    prm = {}
    for line in open(d / "params.txt"):
        k, v = line.split()[:2]
        try:
            prm[k] = float(v)
        except ValueError:
            prm[k] = v
    return prm


def surface_level(s: snap.Snap, r_min_frac: float = 0.7) -> float:
    """外周付近（r > 0.7 Lr）の水柱の高さの平均 = 静かな部分の水面。"""
    i0 = int(s.nr * r_min_frac)
    cols = s.c[i0:, :]
    return float(np.mean(cols.sum(axis=1)) * s.dx)


def analyze_snapshot(s: snap.Snap, prm: dict, jet_on: bool):
    dx = s.dx
    z_s = surface_level(s)
    liq = s.c > 0.5
    st4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    lab, n = ndimage.label(liq, structure=st4)
    R = s.r[:, None] * np.ones((1, s.nz))
    Z = np.ones((s.nr, 1)) * s.z[None, :]
    vol_w = 2 * np.pi * R * s.c * dx * dx
    pool_labels = set(np.unique(lab[:, 0])) - {0}
    top_labels = set(np.unique(lab[:, -1])) - {0} if jet_on else set()
    drops = []
    if n:
        idx = np.arange(1, n + 1)
        vol = ndimage.sum(vol_w, lab, idx)
        zc = ndimage.sum(vol_w * Z, lab, idx) / np.maximum(vol, 1e-30)
        rcm = ndimage.sum(vol_w * R, lab, idx) / np.maximum(vol, 1e-30)
        wc = ndimage.sum(vol_w * s.w, lab, idx) / np.maximum(vol, 1e-30)
        uc = ndimage.sum(vol_w * s.u, lab, idx) / np.maximum(vol, 1e-30)
        cells = ndimage.sum(liq, lab, idx)
        for k in range(n):
            L = k + 1
            if L in pool_labels or L in top_labels:
                continue
            if zc[k] < z_s + 2 * dx:
                continue  # 水中の液体塊（めったにない）
            # 軸上の滴は球、軸外はリング → 3D では周方向に分裂した滴。代表滴径は断面の大きさから
            ring = rcm[k] > 3 * dx
            if ring:
                area = cells[k] * dx * dx
                r_eq = math.sqrt(area / math.pi)  # 断面と同じ面積の円の半径 ≈ 分裂後の滴半径の目安
            else:
                r_eq = (3 * vol[k] / (4 * math.pi)) ** (1 / 3)
            drops.append({"z": float(zc[k]), "r": float(rcm[k]), "vol": float(vol[k]), "w": float(wc[k]),
                          "u": float(uc[k]), "r_eq": float(r_eq), "ring": bool(ring)})
    # プールにつながったまま水面より上に伸びる水（クラウン・柱）
    pool_mask = np.isin(lab, list(pool_labels)) if pool_labels else np.zeros_like(liq)
    a_jet = prm.get("jet_radius", 0.0)
    excl = s.r[:, None] < (a_jet * 1.6 + 2 * dx) if jet_on else np.zeros((s.nr, 1), bool)
    above = pool_mask & (Z > z_s + 2 * dx) & ~excl
    crown_h, crown_w = 0.0, 0.0
    if above.any():
        jj = np.where(above.any(axis=0))[0].max()
        crown_h = (jj + 1) * dx - z_s
        ii = np.where(above[:, jj])[0]
        crown_w = float(np.mean(s.w[ii, jj]))
    # 軸上（水柱がないとき）の盛り上がり = Worthington ジェット
    axis_col = s.c[0, :] > 0.5
    worth_h = 0.0
    if not jet_on:
        js = np.where(axis_col)[0]
        # 水面から連続する部分の最上点
        jz = int(z_s / dx)
        top = jz
        while top + 1 < s.nz and axis_col[top + 1]:
            top += 1
        worth_h = (top + 1) * dx - z_s
    # 巻き込まれた空気（大気につながらない気体の塊）
    air = s.c < 0.5
    alab, an = ndimage.label(air, structure=st4)
    atm = set(np.unique(alab[:, -1])) - {0}
    bub_vol, nb = 0.0, 0
    bubbles = []
    if an:
        aidx = np.arange(1, an + 1)
        avol = ndimage.sum(2 * np.pi * R * (1 - s.c) * dx * dx, alab, aidx)
        azc = ndimage.sum(Z * (1 - s.c), alab, aidx) / np.maximum(ndimage.sum(1 - s.c, alab, aidx), 1e-30)
        arc = ndimage.sum(R * (1 - s.c), alab, aidx) / np.maximum(ndimage.sum(1 - s.c, alab, aidx), 1e-30)
        for k in range(an):
            if (k + 1) in atm:
                continue
            bub_vol += avol[k]; nb += 1
            bubbles.append({"vol": float(avol[k]), "z": float(azc[k]), "r": float(arc[k])})
    return {"t": s.t, "z_s": z_s, "drops": drops, "crown_h": crown_h, "crown_w": crown_w,
            "worth_h": worth_h, "bubble_vol": bub_vol, "n_bubbles": nb, "bubbles": bubbles}


def apex_above_surface(d: dict, z_s: float) -> float:
    """滴の現在位置と速度から、空気抵抗込みで最高到達点（水面基準）を予測。"""
    if d["w"] <= 0:
        return d["z"] - z_s
    sp = math.hypot(d["u"], d["w"])
    ang = math.degrees(math.atan2(d["w"], abs(d["u"]) + 1e-12))
    res = droplets.fly(d["r_eq"], sp, angle_deg=ang)
    return d["z"] - z_s + res["zmax"]


def analyze_run(d: Path, every: int = 1) -> dict:
    prm = load_params(d)
    meta = json.load(open(d / "meta.json")) if (d / "meta.json").exists() else {}
    files = snap.list_snaps(str(d))[::every]
    t_stop = prm.get("jet_stop_time", -1.0)
    rows, all_drops = [], []
    best = {"apex": 0.0}
    for f in files:
        s = snap.read(f)
        jet_on = prm.get("jet_radius", 0) > 0 and (t_stop < 0 or s.t < t_stop)
        a = analyze_snapshot(s, prm, jet_on)
        apexes = [apex_above_surface(dd, a["z_s"]) for dd in a["drops"]]
        rising = [dd for dd in a["drops"] if dd["w"] > 0]
        row = {"t": a["t"], "z_s": a["z_s"], "n_drops": len(a["drops"]), "n_rising": len(rising),
               "drop_vol": sum(dd["vol"] for dd in a["drops"]),
               "max_apex": max(apexes) if apexes else 0.0,
               "max_w": max((dd["w"] for dd in a["drops"]), default=0.0),
               "crown_h": a["crown_h"], "crown_w": a["crown_w"], "worth_h": a["worth_h"],
               "bubble_vol": a["bubble_vol"], "n_bubbles": a["n_bubbles"]}
        rows.append(row)
        for dd, ap in zip(a["drops"], apexes):
            all_drops.append({**dd, "t": a["t"], "apex": ap})
            if ap > best["apex"]:
                best = {"apex": ap, "t": a["t"], **dd}
    # 出力
    import csv
    with open(d / "timeseries.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)
    res = {
        "name": d.name, "meta": meta,
        "max_drop_apex_above_surface": best["apex"], "max_drop": best,
        "max_crown_h": max(r["crown_h"] for r in rows),
        "max_worthington_h": max(r["worth_h"] for r in rows),
        "max_bubble_vol_mm3": max(r["bubble_vol"] for r in rows) * 1e9,
        "final_bubble_vol_mm3": rows[-1]["bubble_vol"] * 1e9,
        "n_drop_detections": len(all_drops),
        "t_end": rows[-1]["t"],
    }
    with open(d / "analysis.json", "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    return res


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = analyze_run(Path(p))
        print(json.dumps({k: v for k, v in r.items() if k != "max_drop"}, indent=1, default=float))
        print("max drop:", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r["max_drop"].items()})
