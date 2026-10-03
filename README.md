# pouring — ウォーターサーバーでコップに水を足すとき、飛沫をどう抑えるか

水が入ったコップに、ウォーターサーバーから水を足す。注いだ瞬間に「水滴の柱」が立ち、
水面では細かい飛沫が跳ねる。ノズルに近づけるとノズルとの距離が近くなりすぎ、
離すと水の勢いがつく。どうするのがいいのか——を、流体シミュレーションで本気で調べたリポジトリです。

## 結論

**コップを傾けて、水流を内壁（水面より 1 cm ほど上）に当てる。水が増えるにつれてコップを起こす。**

1. コップを約 30° 傾け、ノズルの真下に「水面より 1 cm ほど上の内壁」が来るように構える。縁とノズルの隙間は 2〜3 cm。
2. レバーはそっと開ける（0.2 秒くらいかけて）。
3. 水が増えたらコップを起こしていく（335 mL タンブラーの場合: 120 mL で 30°、160 mL で 20°、200 mL で 7°、それ以上は直立）。
4. 直立に戻したら、水流を内壁ぎわに落とす。
5. 傾けられないときは、**水面をノズルから 8 cm 以内に、ただし縁はノズルから 2 cm 以上離す**（10 cm を超えると注ぎ始めの水滴がノズルまで跳ね上がり、泡も急増する。近づけるぶんには泡はほとんどできないので、下限はコップがノズルに触れないことで決まる）。

## 何が起きているのか

| 現象 | 正体 | 効くこと |
|---|---|---|
| 注いだ瞬間の「水滴の柱」 | 水柱の先端の塊が水面に深い穴を掘り、まわりに空気の筒ができる。数十ミリ秒後に筒が潰れる反動で、水柱の外側から細かい水滴が上向きに噴き出す（H = 10 cm で最大約 5 m/s、水面から約 13 cm まで届く計算） | 落下距離を短くする / 壁に当てる / レバーをそっと開ける |
| 水面の細かい飛沫 | 巻き込まれた気泡が浮いてきて弾けるときの「ジェット液滴」。炭酸飲料のミストと同じ。空気抵抗で高さは最大 6 cm 程度 | 空気を巻き込まない（水面に突っ込む速さを、巻き込み開始速度より小さくする） |

空気を巻き込み始める落下距離は、標準的な条件で約 10 cm（水柱の出口での乱れによって 4〜18 cm）。
傾けて壁に当てると、壁が勢いを受け止め、水は滑らかな膜になって斜めに水面へ滑り込むので、桁違いに泡が減ります。

## 3 段構えのシミュレーション

### 1. 軸対称二相流 DNS（`dns/`）
C + OpenMP で書いた、水と空気の非圧縮二相流ソルバー。表面張力まで含めて、注ぎ始めの数百ミリ秒を直接解きます。

- 一様 MAC スタガード格子、軸対称 (r, z)
- 界面: PLIC-VOF（Youngs 法線）+ Weymouth & Yue (2010) の質量保存型の方向分割移流
- 曲率: 高さ関数法（軸対称項込み）、だめなセルは近傍平均 → 平滑化 CSF
- 表面張力: balanced-force CSF
- 圧力: 密度比 ~830 の変係数ポアソン方程式を、Galerkin 幾何マルチグリッド前処理つき Flexible CG で
- 検証: 静止液滴のラプラス圧（誤差 0.1%）、振動液滴の周期（誤差 1.2%）、液滴衝突の Worthington ジェット閾値（Michon ほか 2017 と整合）

```sh
make -C dns                                   # ビルド
python3 dns/cases.py startup --H 0.10         # ノズル〜水面 10 cm の注ぎ始め
python3 dns/batch.py                          # 本番ケース一式
python3 dns/analyze.py dns/runs/startup_H10cm # 水滴の到達点・巻き込み空気を解析
python3 dns/render.py dns/runs/startup_H10cm  # 断面の動画
python3 dns/render3d.py dns/runs/startup_H10cm --t 0.06   # 3D 断面レンダリング
python3 scripts/validate_dns.py               # 検証計算
```

### 2. 半経験物理モデル（`pourphys/`）
文献の式をつないで、あらゆる注ぎ方を一瞬で評価する Python パッケージ。

- `server.py` 吐出（流量・開栓の速さ・ボトルの脈動）
- `jet.py` 落下する水柱の加速・細り・Rayleigh–Plateau 不安定の成長（引き伸ばし込み）
- `startup.py` 注ぎ始めの「先頭の水塊」（粘着粒子モデル = 付着型 Burgers 方程式）
- `bubbles.py` 空気の巻き込み（開始速度・巻き込み量・気泡径）と、気泡破裂のジェット液滴（Deike ほか 2018 の速度則、Gañán-Calvo 型の粒径則）
- `droplets.py` 空気抵抗つきの弾道
- `wall.py` 傾けたコップの内壁を流れる膜
- `cup.py` 傾けたコップの水面・こぼれ限界・水流の当たる点・ノズルとの隙間
- `strategy.py` / `optimize.py` 注ぎ方の評価と探索、水量ごとの傾きの目安

```sh
pip install -r requirements.txt
python3 scripts/sweep_model.py    # results/model.json
python3 -m pytest -q tests
```

### 3. ブラウザで動かせる 2D FLIP（`web/`）
実寸の 2D FLIP 流体（粒子と格子のハイブリッド）でコップと水柱を解き、格子より小さい泡と飛沫は 2 の式で発生させます。
コップをドラッグしたり、傾き・隙間・レバーの開け方・水温を変えたりして試せます。
`web/index.html` をブラウザで開くだけで動きます（レポート付き）。

```sh
python3 scripts/build_web_data.py   # web/data.js と web/media/（DNS の図・動画）を更新
```

## 注意

- 空気を巻き込み始める落下距離は、水柱が出口でどれだけ乱れているかに強く左右されます。
- DNS は軸対称なので、3 次元的な乱れは含みません。軸から離れた「水滴」は計算上はリング状で、実際にはたくさんの粒に分かれます。
- ブラウザ版は 2 次元・表面張力なし。泡と飛沫の数は半経験式による目安です。
- 数値は桁の目安です。傾向（壁に当てると桁違いに減る、近すぎても遠すぎても良くない）は 3 つの方法で一致しています。

## 主な参考文献

- Y. Zhu, H. N. Oğuz, A. Prosperetti, *On the mechanism of air entrainment by liquid jets at a free surface*, J. Fluid Mech. 404 (2000)
- K. T. Kiger, J. H. Duncan, *Air-entrainment mechanisms in plunging jets and breaking waves*, Annu. Rev. Fluid Mech. 44 (2012)
- A. K. Bin, *Gas entrainment by plunging liquid jets*, Chem. Eng. Sci. 48 (1993)
- G.-J. Michon, C. Josserand, T. Séon, *Jet dynamics post drop impact on a deep pool*, Phys. Rev. Fluids 2 (2017)
- L. Deike et al., *Dynamics of jets produced by bursting bubbles*, Phys. Rev. Fluids 3 (2018)
- A. Berny et al., *Role of all jet drops in mass transfer from bursting bubbles*, Phys. Rev. Fluids 5 (2020)
- D. C. Blanchard, *The size and height to which jet drops are ejected from bursting bubbles in seawater*, J. Geophys. Res. 94 (1989)
- G. D. Weymouth, D. K.-P. Yue, *Conservative volume-of-fluid method for free-surface simulations on Cartesian-grid*, J. Comput. Phys. 229 (2010)
- S. Popinet, *An accurate adaptive solver for surface-tension-driven interfacial flows*, J. Comput. Phys. 228 (2009)
