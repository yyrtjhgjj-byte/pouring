"""DNS ケースの生成と実行。

物理的な入力（ノズル〜水面の距離 H、サーバー流量など）から axi2phase のパラメータを作る。

使い方:
  python3 dns/cases.py startup --H 0.10 --res 0.1e-3 --t_end 0.25
  python3 dns/cases.py stop --H 0.10 --t_stop 0.12
  python3 dns/cases.py burst --Rb 0.5e-3
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pourphys.props import Fluid, G  # noqa: E402
from pourphys.server import Server  # noqa: E402

BIN = HERE / "build" / "axi2phase"
RUNS = HERE / "runs"


def fluid_params(fl: Fluid) -> dict:
    return {"rho1": fl.rho, "rho2": fl.rho_air, "mu1": fl.mu, "mu2": fl.mu_air,
            "sigma": fl.sigma, "grav": G}


def grid_for(Lr: float, Lz: float, dx: float, pow2: int = 5) -> tuple[int, int, float]:
    """マルチグリッドが効くよう、2^pow2 の倍数に丸めた格子数。"""
    m = 2 ** pow2
    nr = max(m, int(round(Lr / dx / m)) * m)
    dx = Lr / nr
    nz = max(m, int(round(Lz / dx / m)) * m)
    return nr, nz, dx


def startup_case(H: float, dx: float = 0.1e-3, t_end: float = 0.25, srv: Server | None = None,
                 pool_depth: float = 0.040, L_air: float = 0.030, Lr: float = 0.032,
                 head_factor: float = 1.6, t_stop: float = -1.0, out_dt: float = 1e-3,
                 name: str | None = None) -> dict:
    """注ぎ始め: 先頭の水塊（頭部）付きの水柱が、ノズルから H 下の水面に突っ込む。"""
    srv = srv or Server()
    nr, nz, dx = grid_for(Lr, pool_depth + L_air, dx)
    Lz = nz * dx
    H_top = H - (Lz - pool_depth)          # ノズルからドメイン上端までの落下距離
    if H_top < 0:
        raise ValueError("H が短すぎてドメインに入らない")
    v_top = math.sqrt(srv.v0 ** 2 + 2 * G * H_top)
    a_top = math.sqrt(srv.Q / (math.pi * v_top))
    tip = pool_depth + 0.003
    v_tip = math.sqrt(v_top ** 2 + 2 * G * (Lz - tip))
    a_tip = a_top * math.sqrt(v_top / v_tip)
    prm = {
        "nr": nr, "nz": nz, "Lr": Lr, **fluid_params(srv.fluid),
        "pool_depth": pool_depth, "jet_radius": a_top, "jet_velocity": v_top,
        "jet_head_radius": head_factor * a_tip, "jet_tip_z": tip, "jet_init": 1, "jet_profile": 1,
        "jet_stop_time": t_stop, "t_end": t_end, "out_dt": out_dt, "tol": 1e-4, "cfl": 0.4,
    }
    name = name or (f"startup_H{H*100:.0f}cm" + (f"_stop{t_stop*1000:.0f}ms" if t_stop > 0 else ""))
    meta = {"kind": "startup", "H": H, "U_impact": v_tip, "a_impact": a_tip, "v_top": v_top,
            "a_top": a_top, "Lz": Lz, "dx": dx, "head_D": 2 * head_factor * a_tip}
    return {"name": name, "params": prm, "meta": meta}


def drop_case(D: float, Fr: float, dx: float = 0.0625e-3, t_end: float = 0.08,
              fl: Fluid | None = None, name: str | None = None) -> dict:
    fl = fl or Fluid.water(20.0)
    U = math.sqrt(Fr * G * D)
    Lr, pool, air = 6 * D, 6 * D, 8 * D
    nr, nz, dx = grid_for(Lr, pool + air, dx)
    prm = {"nr": nr, "nz": nz, "Lr": Lr, **fluid_params(fl), "pool_depth": pool,
           "drop_radius": D / 2, "drop_z": pool + D / 2 + 3 * dx, "drop_velocity": U,
           "t_end": t_end, "out_dt": 0.5e-3, "tol": 1e-4}
    return {"name": name or f"drop_D{D*1e3:.1f}mm_Fr{Fr:.0f}", "params": prm,
            "meta": {"kind": "drop", "D": D, "Fr": Fr, "U": U, "pool": pool}}


def burst_case(Rb: float, cells_per_R: int = 48, t_end: float | None = None,
               fl: Fluid | None = None, name: str | None = None) -> dict:
    """気泡破裂: 水面直下の球状のくぼみ（膜はすでに破れた状態）から始める。"""
    fl = fl or Fluid.water(20.0)
    Lr, pool, air = 4 * Rb, 4 * Rb, 10 * Rb
    nr, nz, dx = grid_for(Lr, pool + air, Rb / cells_per_R)
    open_r = 0.25 * Rb
    depth = math.sqrt(Rb ** 2 - open_r ** 2)
    tc = math.sqrt(fl.rho * Rb ** 3 / fl.sigma)
    prm = {"nr": nr, "nz": nz, "Lr": Lr, **fluid_params(fl), "pool_depth": pool,
           "cavity_radius": Rb, "cavity_depth": depth, "t_end": t_end or 4 * tc,
           "out_dt": tc / 40, "tol": 1e-4, "cfl": 0.3}
    return {"name": name or f"burst_R{Rb*1e6:.0f}um", "params": prm,
            "meta": {"kind": "burst", "Rb": Rb, "pool": pool, "t_cap": tc}}


def write_params(case: dict, root: Path = RUNS) -> Path:
    d = root / case["name"]
    d.mkdir(parents=True, exist_ok=True)
    prm = dict(case["params"])
    prm["out_dir"] = str(d)
    with open(d / "params.txt", "w") as f:
        for k, v in prm.items():
            f.write(f"{k} {v}\n")
    import json
    with open(d / "meta.json", "w") as f:
        json.dump(case["meta"], f, indent=1)
    return d


def run(case: dict, root: Path = RUNS, threads: int | None = None) -> Path:
    if not BIN.exists():
        subprocess.check_call(["make", "-C", str(HERE)])
    d = write_params(case, root)
    env = dict(os.environ)
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
    with open(d / "log.txt", "w") as log:
        subprocess.check_call([str(BIN), str(d / "params.txt")], stderr=log, env=env)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["startup", "stop", "drop", "burst"])
    ap.add_argument("--H", type=float, default=0.10)
    ap.add_argument("--res", type=float, default=0.1e-3)
    ap.add_argument("--t_end", type=float, default=None)
    ap.add_argument("--t_stop", type=float, default=0.12)
    ap.add_argument("--D", type=float, default=2.6e-3)
    ap.add_argument("--Fr", type=float, default=100.0)
    ap.add_argument("--Rb", type=float, default=0.5e-3)
    ap.add_argument("--name", default=None)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    if a.kind == "startup":
        case = startup_case(a.H, dx=a.res, t_end=a.t_end or 0.25, name=a.name)
    elif a.kind == "stop":
        case = startup_case(a.H, dx=a.res, t_end=a.t_end or (a.t_stop + 0.15), t_stop=a.t_stop, name=a.name)
    elif a.kind == "drop":
        case = drop_case(a.D, a.Fr, t_end=a.t_end or 0.08, name=a.name)
    else:
        case = burst_case(a.Rb, t_end=a.t_end, name=a.name)
    print(case["name"], case["meta"])
    if a.dry:
        print(write_params(case))
    else:
        print("->", run(case))


if __name__ == "__main__":
    main()
