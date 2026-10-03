"""DNS ソルバーの検証。結果は results/dns_validation.json。

1. 静止液滴: ラプラス圧 Δp = 2σ/R と寄生流
2. 振動液滴: n=2 モードの周期（Lamb の理論）
3. 液滴の深いプールへの衝突: Worthington ジェットが出る Fr の閾値（Michon et al. 2017: Fr ≈ 60 / 90）

  python3 scripts/validate_dns.py
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "dns"))
import snap  # noqa: E402

BIN = ROOT / "dns" / "build" / "axi2phase"
VDIR = ROOT / "dns" / "runs" / "validation"
SIG, RHO, RHOG = 0.0727, 998.2, 1.204


def run(name, prm):
    d = VDIR / name
    d.mkdir(parents=True, exist_ok=True)
    prm = dict(prm, out_dir=str(d), verbose=0)
    with open(d / "params.txt", "w") as f:
        for k, v in prm.items():
            f.write(f"{k} {v}\n")
    if not ((d / "log.txt").exists() and "done:" in (d / "log.txt").read_text()):
        with open(d / "log.txt", "w") as log:
            subprocess.check_call([str(BIN), str(d / "params.txt")], stderr=log)
    return d


def main():
    out = {}
    # 1. 静止液滴
    R = 0.002
    d = run("static_drop", {"nr": 64, "nz": 128, "Lr": 0.004, "grav": 0, "pool_depth": 0, "sigma": SIG,
                            "drop_radius": R, "drop_z": 0.004, "t_end": 0.01, "out_dt": 0.002})
    s = snap.read(snap.list_snaps(str(d))[-1])
    dp = float(s.p[s.c > 0.999].mean() - s.p[s.c < 0.001].mean())
    out["static_drop"] = {"R_mm": R * 1e3, "cells_per_R": 32, "dp_sim": dp, "dp_theory": 2 * SIG / R,
                          "rel_err": dp / (2 * SIG / R) - 1, "max_spurious_velocity": float(np.hypot(s.u, s.w).max()),
                          "capillary_number": float(1e-3 * np.hypot(s.u, s.w).max() / SIG)}
    # 2. 振動液滴
    d = run("osc_drop", {"nr": 64, "nz": 128, "Lr": 0.006, "grav": 0, "pool_depth": 0, "sigma": SIG,
                         "drop_radius": R, "drop_z": 0.006, "drop_p2": 0.08, "t_end": 0.05, "out_dt": 0.0005})
    ts, hs = [], []
    for f in snap.list_snaps(str(d)):
        s = snap.read(f)
        ts.append(s.t); hs.append(s.c[0].sum() * s.dx / 2)
    ts, hs = np.array(ts), np.array(hs) - np.mean(hs)
    zc = [ts[i] - hs[i] * (ts[i + 1] - ts[i]) / (hs[i + 1] - hs[i]) for i in range(len(hs) - 1) if hs[i] * hs[i + 1] < 0]
    per = 2 * float(np.mean(np.diff(zc)))
    w2 = 24 * SIG / (R ** 3 * (3 * RHO + 2 * RHOG))
    out["oscillating_drop"] = {"period_sim_ms": per * 1e3, "period_theory_ms": 2 * math.pi / math.sqrt(w2) * 1e3,
                               "rel_err": per / (2 * math.pi / math.sqrt(w2)) - 1}
    # 3. 液滴衝突
    D = 0.0026
    rows = []
    for Fr in (40, 70, 110, 160):
        U = math.sqrt(Fr * 9.80665 * D)
        d = run(f"impact_Fr{Fr}", {"nr": 192, "nz": 384, "Lr": 0.012, "pool_depth": 0.014, "sigma": SIG,
                                   "drop_radius": D / 2, "drop_z": 0.0155, "drop_velocity": U,
                                   "t_end": 0.08, "out_dt": 0.001, "tol": 1e-4})
        hmax, tmax = 0.0, 0.0
        for f in snap.list_snaps(str(d)):
            s = snap.read(f)
            if s.t < 0.004:
                continue
            j = np.where(s.c[0] > 0.5)[0]
            h = (j.max() + 1) * s.dx - 0.014 if len(j) else 0
            if h > hmax:
                hmax, tmax = h, s.t
        rows.append({"Fr": Fr, "U": U, "jet_height_mm": hmax * 1e3, "t_ms": tmax * 1e3,
                     "regime": "jet" if hmax > 0.004 else "no jet"})
    out["drop_impact"] = {"D_mm": D * 1e3, "rows": rows,
                          "reference": "Michon, Josserand & Séon, Phys. Rev. Fluids 2, 023601 (2017): Fr≲60 ジェットなし, Fr≳90 高速ジェット"}
    (ROOT / "results").mkdir(exist_ok=True)
    with open(ROOT / "results" / "dns_validation.json", "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
