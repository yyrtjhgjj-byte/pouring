"""レポート用のデータと図をまとめる。

  python3 scripts/build_web_data.py            # web/data.js と web/media/*
  python3 scripts/build_web_data.py --no-media # data.js だけ

- results/model.json（半経験モデルの掃引）
- dns/runs/*/analysis.json（DNS の解析; なければ dns/analyze.py を実行）
- DNS のモンタージュ画像と動画
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "dns"))
sys.path.insert(0, str(ROOT))
RUNS = ROOT / "dns" / "runs"
WEB = ROOT / "web"
MEDIA = WEB / "media"


def run_done(d: Path) -> bool:
    log = d / "log.txt"
    return log.exists() and "done:" in log.read_text()


def ensure_analysis(d: Path) -> dict | None:
    if not run_done(d):
        return None
    a = d / "analysis.json"
    if not a.exists() or a.stat().st_mtime < (d / "log.txt").stat().st_mtime:
        import analyze
        analyze.analyze_run(d)
    return json.load(open(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-media", action="store_true")
    a = ap.parse_args()
    model = json.load(open(ROOT / "results" / "model.json"))
    dns = {"startup": [], "stop": [], "burst": [], "validation": {}}
    for d in sorted(RUNS.glob("startup_H*cm*")):
        res = ensure_analysis(d)
        if not res:
            continue
        meta = res["meta"]
        row = {"name": d.name, "H": meta["H"], "U": meta["U_impact"], "head_D_mm": meta["head_D"] * 1e3,
               "apex": res["max_drop_apex_above_surface"], "air_mL": res["max_bubble_vol_mm3"] / 1000,
               "air_final_mL": res["final_bubble_vol_mm3"] / 1000, "crown_mm": res["max_crown_h"] * 1e3,
               "worth_mm": res["max_worthington_h"] * 1e3, "t_apex": res["max_drop"].get("t", None),
               "w_max": res["max_drop"].get("w", 0.0), "drop_r_um": res["max_drop"].get("r_eq", 0) * 1e6,
               "apex_p90": res.get("apex_p90", 0.0), "apex_p50": res.get("apex_p50", 0.0),
               "n_above_nozzle": res.get("n_above_nozzle", 0),
               "stop_apex": res.get("stop_apex", 0.0), "stop_worth_mm": res.get("max_worthington_after_stop", 0.0) * 1e3}
        if "_stop" in d.name:
            # 止めた直後は水柱の尻尾が軸上に残るので、尻尾が吸い込まれた後（停止 50 ms 後以降）の盛り上がりだけ数える
            import csv
            t_stop = float([l.split()[1] for l in open(d / "params.txt") if l.startswith("jet_stop_time")][0])
            ts = list(csv.DictReader(open(d / "timeseries.csv")))
            row["t_stop"] = t_stop
            row["worth_mm"] = max((float(r["worth_h"]) for r in ts if float(r["t"]) > t_stop + 0.05), default=0.0) * 1e3
            dns["stop"].append(row)
        elif "nohead" in d.name:
            dns["nohead"] = row
        else:
            dns["startup"].append(row)
    from pourphys import bubbles as pb, ROOM
    for d in sorted(RUNS.glob("burst_R*")):
        res = ensure_analysis(d)
        if not res or not res.get("first_axis_drop"):
            continue
        Rb = res["meta"]["Rb"]
        rd, Vd, n = pb.jet_drops(Rb, ROOM)
        fd = res["first_axis_drop"]
        dns["burst"].append({"name": d.name, "Rb": Rb, "v_sim": fd["w"], "r_sim_um": fd["r_eq"] * 1e6,
                             "t_ms": fd["t"] * 1e3, "v_theory": float(Vd), "r_theory_um": float(rd) * 1e6})
    # 図（存在するものだけ）
    frames = [("r3d_stop14_160ms.png", "160 ms · 尻尾が吸い込まれる", "流れを止めて 10 ミリ秒後。水柱の尻尾が円錐状に細くなって水面に吸い込まれる"),
              ("r3d_stop14_190ms.png", "190 ms · 真ん中に小さな柱", "真ん中に小さな柱が立ち、先端から粒がちぎれる"),
              ("r3d_stop14_200ms.png", "200 ms · 粒が真上へ", "ちぎれた粒が真上に飛ぶ"),
              ("r3d_stop14_215ms.png", "215 ms · 水滴の列", "真ん中の上に水滴が一列に並んで上がっていく")]
    dns["stop_frames"] = [{"src": f"media/{f}", "cap": c, "alt": al} for f, c, al in frames if (MEDIA / f).exists()]
    bframes = sorted(MEDIA.glob("burst_R*_*.png"))
    dns["burst_frames"] = [{"src": f"media/{f.name}", "cap": f.stem.split("_")[-1].replace("us", " µs").replace("ms", " ms"),
                            "alt": "泡の破裂の断面"} for f in bframes]
    vfile = ROOT / "results" / "dns_validation.json"
    if vfile.exists():
        dns["validation"] = json.load(open(vfile))
    data = {"model": model, "dns": dns}
    WEB.mkdir(exist_ok=True)
    with open(WEB / "data.js", "w") as f:
        f.write("// 自動生成: scripts/build_web_data.py\nwindow.POUR_DATA = ")
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"), default=float)
        f.write(";\n")
    print("wrote web/data.js;", len(dns["startup"]), "startup runs,", len(dns["stop"]), "stop runs,", len(dns["burst"]), "burst runs")
    if a.no_media:
        return
    MEDIA.mkdir(exist_ok=True)
    import render
    for row in dns["startup"] + dns["stop"] + ([dns["nohead"]] if "nohead" in dns else []):
        d = RUNS / row["name"]
        mp4 = d / "movie.mp4"
        if not mp4.exists():
            render.render_movie(d, every=1, zoom=(-25, 25, 12, d and 70.4), dpi=90)
        shutil.copy(mp4, MEDIA / f"{row['name']}.mp4")
    print("media copied")


if __name__ == "__main__":
    main()
