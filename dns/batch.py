"""本番 DNS をまとめて順番に実行する。

  python3 dns/batch.py            # 全ケース
  python3 dns/batch.py --only startup_H10cm
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cases  # noqa: E402


def all_cases(dx: float = 0.1e-3):
    cs = []
    # 1) 注ぎ始め: ノズル〜水面の距離 H を変える（先頭の水塊 = 水柱径の 1.6 倍）
    for H in [0.10, 0.04, 0.07, 0.14, 0.20]:
        cs.append(cases.startup_case(H, dx=dx, t_end=0.20))
    # 2) 先頭の水塊の大きさの感度（丸い先端だけ = 水柱径と同じ）
    cs.append(cases.startup_case(0.10, dx=dx, t_end=0.20, head_factor=1.0, name="startup_H10cm_nohead"))
    # 3) 注ぎ終わり（流入停止）→ 尾部の Worthington ジェット
    for H in [0.07, 0.14]:
        cs.append(cases.startup_case(H, dx=dx, t_end=0.30, t_stop=0.15))
    # 4) 気泡破裂（微細飛沫の検証）
    for Rb in [0.5e-3, 0.25e-3]:
        cs.append(cases.burst_case(Rb))
    return cs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    for case in all_cases():
        if a.only and case["name"] not in a.only:
            continue
        d = cases.RUNS / case["name"]
        if a.skip_done and (d / "log.txt").exists() and "done:" in (d / "log.txt").read_text():
            print("skip", case["name"], flush=True)
            continue
        t0 = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] start {case['name']}", flush=True)
        cases.run(case)
        print(f"[{time.strftime('%H:%M:%S')}] done  {case['name']} in {time.time()-t0:.0f}s", flush=True)
