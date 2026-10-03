/*
 * sim.js — ウォーターサーバーの注水シミュレーター（ブラウザ／Node 両対応のコア）
 *
 *  ・2D FLIP 流体（粒子＋格子のハイブリッド, 実寸 [m] で計算）
 *      水柱・コップの水・コップ（動かせる剛体の壁）を解く
 *  ・サブグリッド物理（格子より小さい現象は pourphys と同じ半経験式で）
 *      水柱が水面／壁に当たった瞬間の速さ・表面の乱れ → 空気の巻き込み → 気泡
 *      → 浮上して水面で破裂 → ジェット液滴（微細飛沫）→ 空気抵抗つき弾道
 *  ・ノズルに当たった飛沫、縁より上に飛んだ飛沫を数える
 *
 * 座標: x 右, y 上, 単位 m。原点はドメイン左下。
 */
(function (root) {
  'use strict';
  const G = 9.80665;

  // ======================================================================
  // 物性・半経験式（pourphys の移植）
  // ======================================================================
  function waterProps(T) {
    const a1 = -3.983035, a2 = 301.797, a3 = 522528.9, a4 = 69.34881, a5 = 999.97495;
    const rho = a5 * (1 - (T + a1) * (T + a1) * (T + a2) / (a3 * (T + a4)));
    const TK = T + 273.15;
    const mu = 2.414e-5 * Math.pow(10, 247.8 / (TK - 140));
    const tau = 1 - TK / 647.096;
    const sigma = 235.8e-3 * Math.pow(tau, 1.256) * (1 - 0.625 * tau);
    return { T, rho, mu, sigma, nu: mu / rho, lmu: mu * mu / (rho * sigma) };
  }

  // 第 1 種変形ベッセル関数の比 I1(x)/I0(x)（級数, x<=1.5 で十分）
  function besselRatio(x) {
    let i0 = 1, i1 = x / 2, t0 = 1, t1 = x / 2;
    for (let k = 1; k < 12; k++) {
      t0 *= (x * x / 4) / (k * k); i0 += t0;
      t1 *= (x * x / 4) / (k * (k + 1)); i1 += t1;
    }
    return i1 / i0;
  }

  // 落下ジェットの表面乱れの成長（Rayleigh–Plateau, 引き伸ばし込み）
  function jetTable(Q, d0, fl, zmax = 0.4, n = 240, nk = 60) {
    const a0 = d0 / 2, v0 = Q / (Math.PI * a0 * a0);
    const z = new Float64Array(n), v = new Float64Array(n), a = new Float64Array(n), t = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      z[i] = zmax * i / (n - 1);
      v[i] = Math.sqrt(v0 * v0 + 2 * G * z[i]);
      a[i] = Math.sqrt(Q / (Math.PI * v[i]));
      if (i > 0) t[i] = t[i - 1] + (z[i] - z[i - 1]) * 0.5 * (1 / v[i] + 1 / v[i - 1]);
    }
    const gainBest = new Float64Array(n);
    const gk = new Float64Array(nk);
    let prevOm = new Float64Array(nk);
    for (let i = 0; i < n; i++) {
      let best = 0;
      for (let k = 0; k < nk; k++) {
        const x0 = 0.05 + 1.45 * k / (nk - 1);
        const x = (x0 / a0) * (v0 / v[i]) * a[i];
        const cap = fl.sigma / (fl.rho * a[i] * a[i] * a[i]);
        const wR2 = x < 1 ? cap * x * (1 - x * x) * besselRatio(x) : 0;
        const s = 1.5 * fl.nu * x * x / (a[i] * a[i]);
        const om = Math.sqrt(wR2 + s * s) - s;
        if (i > 0) gk[k] += 0.5 * (om + prevOm[k]) * (t[i] - t[i - 1]);
        prevOm[k] = om;
        if (gk[k] > best) best = gk[k];
      }
      gainBest[i] = best;
    }
    return { z, v, a, gain: gainBest, v0, zmax };
  }
  function interpTable(tab, arr, zq) {
    const n = tab.z.length;
    const f = Math.min(Math.max(zq / tab.zmax, 0), 1) * (n - 1);
    const i = Math.min(Math.floor(f), n - 2), w = f - i;
    return arr[i] * (1 - w) + arr[i + 1] * w;
  }

  const V_E_SMOOTH = 3.0, V_E_ROUGH = 0.8, DELTA_REF = 0.08, LA_C = 500;
  function onsetVelocity(rough, fl) {
    return Math.sqrt(fl.sigma / 0.0728) * (V_E_ROUGH + (V_E_SMOOTH - V_E_ROUGH) * Math.exp(-rough / DELTA_REF));
  }
  function airRatio(U, H, d0, v0, rough, fl) {
    const Ve = onsetVelocity(rough, fl);
    const Fr0 = v0 * v0 / (G * d0);
    const base = 0.04 * Math.pow(Fr0, 0.28) * Math.pow(Math.max(H / d0, 1), 0.4);
    const x = U / Ve;
    const gate = Math.pow(Math.max(1 - 1 / (x * x), 0), 1.5) + 0.02 / (1 + Math.exp(-(x - 1) / 0.06));
    return base * Math.min(gate, 1);
  }
  function bubbleRmed(U, d, fl) {
    const eps = 0.05 * U * U * U / Math.max(d, 1e-4);
    const dH = 0.725 * Math.pow(fl.sigma / fl.rho, 0.6) * Math.pow(Math.max(eps, 1e-6), -0.4);
    return Math.min(Math.max(0.5 * dH, 30e-6), 2e-3);
  }
  function jetDrops(R, fl) {
    const La = R / fl.lmu;
    const Bo = fl.rho * G * R * R / fl.sigma;
    if (La < LA_C) return { rd: 0, Vd: 0, n: 0 };
    const Vc = Math.sqrt(fl.sigma / (fl.rho * R));
    const s = Math.sqrt(La / LA_C) - 1 + 1e-9;
    const Vd = 19 * Vc * Math.pow(s, -0.5) * Math.pow(1 + 2.2 * Bo, -0.75);
    const rd = R * 0.6 * Math.pow(s, 1.25) / Math.pow(La, 0.375);
    const n = Math.max(0, Math.min(7, Math.round(7 * (1 - Bo / 0.3))));
    return { rd, Vd, n };
  }
  function riseVelocity(R, fl) {
    // 小さい気泡は（汚れた水道水想定で）剛体球 Stokes、上限 0.22 m/s
    return Math.min(2 / 9 * G * R * R / fl.nu, 0.22);
  }
  // 正規乱数
  function randn() {
    let u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  // ======================================================================
  // コップ（剛体）: 3 本のカプセル（左壁・右壁・底）
  // ======================================================================
  class Cup {
    constructor(o = {}) {
      this.rb = o.rb ?? 0.030;   // 底の内半幅
      this.rt = o.rt ?? 0.037;   // 縁の内半幅
      this.hc = o.hc ?? 0.095;   // 内高さ
      this.tw = o.tw ?? 0.0034;  // 壁の厚み
      this.x = o.x ?? 0.08; this.y = o.y ?? 0.02; this.th = o.th ?? 0; // 内底中心の位置と傾き
      this.vx = 0; this.vy = 0; this.om = 0;
      this.update();
    }
    // ローカル → ワールド
    toWorld(lx, ly) {
      const c = Math.cos(this.th), s = Math.sin(this.th);
      return [this.x + lx * c - ly * s, this.y + lx * s + ly * c];
    }
    toLocal(wx, wy) {
      const c = Math.cos(this.th), s = Math.sin(this.th);
      const dx = wx - this.x, dy = wy - this.y;
      return [dx * c + dy * s, -dx * s + dy * c];
    }
    update() {
      const r = this.tw / 2, rb = this.rb + r, rt = this.rt + r;
      const segL = [-rb, -r, -rt, this.hc], segR = [rb, -r, rt, this.hc], segB = [-rb, -r, rb, -r];
      this.segs = [segL, segR, segB].map(([ax, ay, bx, by]) => {
        const A = this.toWorld(ax, ay), B = this.toWorld(bx, by);
        return [A[0], A[1], B[0], B[1]];
      });
      this.rad = r;
      let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
      for (const sg of this.segs) {
        x0 = Math.min(x0, sg[0], sg[2]); x1 = Math.max(x1, sg[0], sg[2]);
        y0 = Math.min(y0, sg[1], sg[3]); y1 = Math.max(y1, sg[1], sg[3]);
      }
      this.bbox = [x0 - 0.01, x1 + 0.01, y0 - 0.01, y1 + 0.01];
      // 縁の点
      this.rimL = this.toWorld(-this.rt - r, this.hc);
      this.rimR = this.toWorld(this.rt + r, this.hc);
    }
    wallVel(px, py) {
      return [this.vx - this.om * (py - this.y), this.vy + this.om * (px - this.x)];
    }
    // 割り当てなしの高速版: this._d, this._nx, this._ny に結果を書く
    dist2(px, py) {
      let best = 1e9, bx = 0, by = 1;
      const S = this.segs;
      for (let q = 0; q < 3; q++) {
        const sg = S[q];
        const ax = sg[0], ay = sg[1], ex = sg[2] - ax, ey = sg[3] - ay;
        let t = ((px - ax) * ex + (py - ay) * ey) / (ex * ex + ey * ey);
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const dx = px - ax - t * ex, dy = py - ay - t * ey, d2 = dx * dx + dy * dy;
        if (d2 < best) { best = d2; bx = dx; by = dy; }
      }
      const d = Math.sqrt(best);
      this._d = d; this._nx = d > 1e-12 ? bx / d : 0; this._ny = d > 1e-12 ? by / d : 1;
      return d;
    }
    // 点からの最短距離（壁の芯線まで）と法線
    distance(px, py) {
      let best = 1e9, nx = 0, ny = 0;
      for (const [ax, ay, bx, by] of this.segs) {
        const ex = bx - ax, ey = by - ay;
        let t = ((px - ax) * ex + (py - ay) * ey) / (ex * ex + ey * ey);
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const qx = ax + t * ex, qy = ay + t * ey;
        const dx = px - qx, dy = py - qy, d = Math.hypot(dx, dy);
        if (d < best) { best = d; nx = d > 1e-12 ? dx / d : 0; ny = d > 1e-12 ? dy / d : 1; }
      }
      return [best, nx, ny];
    }
    inside(px, py) {
      // 内側（水が入る空間）か
      const [lx, ly] = this.toLocal(px, py);
      if (ly < 0 || ly > this.hc) return false;
      const hw = this.rb + (this.rt - this.rb) * ly / this.hc;
      return Math.abs(lx) < hw;
    }
    capacityArea() { return (this.rb + this.rt) * this.hc; }
  }

  // ======================================================================
  // FLIP シミュレーター
  // ======================================================================
  const FLUID = 0, AIR = 1, SOLID = 2;

  class PourSim {
    constructor(cfg = {}) {
      this.cfg = cfg;
      this.h = cfg.h ?? 0.00125;
      this.nx = Math.floor((cfg.W ?? 0.26) / this.h);
      this.ny = Math.floor((cfg.H ?? 0.28) / this.h);
      this.W = this.nx * this.h; this.H = this.ny * this.h;
      const n = this.nx * this.ny;
      this.u = new Float32Array(n); this.v = new Float32Array(n);
      this.du = new Float32Array(n); this.dv = new Float32Array(n);
      this.pu = new Float32Array(n); this.pv = new Float32Array(n);
      this.s = new Float32Array(n); this.type = new Int32Array(n);
      this.p = new Float32Array(n);
      this.dens = new Float32Array(n); this.poolCnt = new Float32Array(n); this.filmCnt = new Float32Array(n);
      this.solidVx = new Float32Array(n); this.solidVy = new Float32Array(n);
      this.restDensity = 0;
      // 粒子
      this.pr = 0.3 * this.h;                 // 粒子半径
      this.dp = this.h / 2;                   // 初期粒子間隔（1 セル 4 粒子）
      this.cap = cfg.maxParticles ?? 40000;
      this.pos = new Float32Array(2 * this.cap); this.vel = new Float32Array(2 * this.cap);
      this.state = new Int8Array(this.cap);   // 0: 水柱(空中) 1: 壁膜 2: 水だまり
      this.yb = new Float32Array(this.cap);   // 放出時の y（落下距離用）
      this.rough = new Float32Array(this.cap);
      this.yWall = new Float32Array(this.cap); // 壁に触れた高さ
      this.np = 0;
      // 近傍探索
      this.pInv = 1 / (2.2 * this.pr);
      this.pnx = Math.floor(this.W * this.pInv) + 1; this.pny = Math.floor(this.H * this.pInv) + 1;
      this.cellCount = new Int32Array(this.pnx * this.pny + 1);
      this.cellIds = new Int32Array(this.cap);
      // サーバー
      this.xN = cfg.xN ?? this.W * 0.40;
      this.yN = cfg.yN ?? this.H - 0.022;
      this.d0 = cfg.d0 ?? 0.007;
      this.rNoz = 0.0075;
      this.Q = cfg.Q ?? 25e-6;
      this.tauValve = cfg.tauValve ?? 0.06;
      this.T = cfg.T ?? 20;
      this.fl = waterProps(this.T);
      this.noise = 0.02;
      this.rebuildTables();
      this.lever = false; this.leverT = 0; this.flowFrac = 0; this.emitAcc = 0;
      this.levelHold = cfg.levelHold ?? true;   // 注いだ分だけコップの底から抜いて水位を保つ（2D は実物より速く満杯になるため）
      this.emittedSinceDrain = 0;
      // コップ
      this.cup = new Cup(cfg.cup || {});
      this.cup.x = cfg.cupX ?? this.xN; this.cup.y = cfg.cupY ?? 0.02;
      this.cup.update();
      // サブグリッド: 気泡・飛沫
      this.maxB = 2500; this.maxD = 3000;
      this.bx = new Float32Array(this.maxB); this.by = new Float32Array(this.maxB);
      this.bR = new Float32Array(this.maxB); this.bw = new Float32Array(this.maxB); this.nb = 0;
      this.dx_ = new Float32Array(this.maxD); this.dy_ = new Float32Array(this.maxD);
      this.dvx = new Float32Array(this.maxD); this.dvy = new Float32Array(this.maxD);
      this.dr = new Float32Array(this.maxD); this.dw = new Float32Array(this.maxD);
      this.dTop = new Float32Array(this.maxD); this.dFlag = new Uint8Array(this.maxD); this.nd = 0;
      this.time = 0;
      this.resetStats();
      this.flip = 0.95;
      this.pressureIters = cfg.pressureIters ?? 30;
      this.markSolids();
      if (cfg.fill !== false) this.fillCup(cfg.fillFrac ?? 0.38);
    }

    rebuildTables() {
      this.fl = waterProps(this.T);
      this.tab = jetTable(this.Q, this.d0, this.fl);
      this.v0 = this.Q / (Math.PI * this.d0 * this.d0 / 4);
      // 3D 相当の 1 粒子あたりの体積: 定常放出率 = v0 d0 / dp^2 [個/s]
      this.Vp = this.Q * this.dp * this.dp / (this.v0 * this.d0);
    }
    setTemperature(T) { this.T = T; this.rebuildTables(); }
    setFlow(Q) { this.Q = Q; this.rebuildTables(); }

    resetStats() {
      this.stats = {
        nozzleHits: 0, rimEscape: 0, bigNozzleHits: 0, spilled: 0,
        bubbles: 0, drops: 0,
        // 直近の着水情報（指数移動平均）
        H: 0, U: 0, Ueff: 0, rough: 0, Ve: 0, air: 0, wallFrac: 0, hasImpact: false,
        bubbleRate: 0, dropRate: 0, events: [],
      };
      this._accB = 0; this._accD = 0; this._accT = 0;
    }

    idx(i, j) { return i * this.ny + j; }

    // ---------------------------------------------------------- 固体セル
    markSolids() {
      const { nx, ny, h, cup } = this;
      cup.update();
      const R = cup.rad + 0.35 * h;
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        const k = i * ny + j;
        let solid = (i === 0 || i === nx - 1 || j === 0);
        const x = (i + 0.5) * h, y = (j + 0.5) * h;
        if (!solid) {
          const [d] = cup.distance(x, y);
          if (d < R) {
            solid = true;
            const [wx, wy] = cup.wallVel(x, y);
            this.solidVx[k] = wx; this.solidVy[k] = wy;
          }
        }
        // サーバーの本体（ノズルより上、出口を除く）
        if (!solid && y > this.yN + 0.004 && Math.abs(x - this.xN) < 0.03) solid = true;
        if (!solid && y > this.yN && Math.abs(x - this.xN) < this.rNoz && Math.abs(x - this.xN) > this.d0 / 2) solid = true;
        if (solid && !(cup.distance(x, y)[0] < R)) { this.solidVx[k] = 0; this.solidVy[k] = 0; }
        this.s[k] = solid ? 0 : 1;
      }
    }

    fillCup(frac) {
      const { cup, dp } = this;
      const yTop = cup.hc * frac;
      for (let ly = dp / 2; ly < yTop; ly += dp)
        for (let lx = -cup.rt; lx <= cup.rt; lx += dp) {
          const hw = cup.rb + (cup.rt - cup.rb) * ly / cup.hc - this.pr * 1.2;
          if (Math.abs(lx) > hw) continue;
          const [x, y] = cup.toWorld(lx + (Math.random() - 0.5) * 1e-5, ly);
          this.addParticle(x, y, 0, 0, 2, y);
        }
      this.transfer(true);
      this.updateDensity();
      // 静止密度（水の中の平均）
      let sum = 0, cnt = 0;
      for (let k = 0; k < this.dens.length; k++) if (this.dens[k] > 0 && this.type[k] === FLUID) { sum += this.dens[k]; cnt++; }
      this.restDensity = cnt ? sum / cnt : 4;
    }

    addParticle(x, y, vx, vy, st, yb) {
      if (this.np >= this.cap) return -1;
      const i = this.np++;
      this.pos[2 * i] = x; this.pos[2 * i + 1] = y; this.vel[2 * i] = vx; this.vel[2 * i + 1] = vy;
      this.state[i] = st; this.yb[i] = yb; this.rough[i] = 0;
      return i;
    }
    removeParticle(i) {
      const j = --this.np;
      if (i === j) return;
      this.pos[2 * i] = this.pos[2 * j]; this.pos[2 * i + 1] = this.pos[2 * j + 1];
      this.vel[2 * i] = this.vel[2 * j]; this.vel[2 * i + 1] = this.vel[2 * j + 1];
      this.state[i] = this.state[j]; this.yb[i] = this.yb[j]; this.rough[i] = this.rough[j]; this.yWall[i] = this.yWall[j];
    }

    // ---------------------------------------------------------- 吐出
    exitVelocity() {
      return this.v0 * this.flowFrac;
    }
    emit(dt) {
      if (this.lever) this.leverT += dt; else this.leverT = 0;
      const target = this.lever ? 1 : 0;
      // 開栓は smoothstep で τ_v かけて立ち上がる。閉栓は速い（0.03 s）
      if (this.lever) {
        const s = Math.min(this.leverT / Math.max(this.tauValve, 1e-3), 1);
        this.flowFrac = s * s * (3 - 2 * s);
      } else {
        this.flowFrac = Math.max(0, this.flowFrac - dt / 0.03);
      }
      const v = this.exitVelocity();
      if (v <= 1e-4) return;
      this.emitAcc += v * dt;
      const dp = this.dp;
      const nAcross = Math.max(1, Math.round(this.d0 / dp));
      // 低流量のときは細く出る（満管でない）
      const width = this.d0 * Math.sqrt(Math.max(this.flowFrac, 0.05));
      while (this.emitAcc >= dp) {
        this.emitAcc -= dp;
        const yy = this.yN - this.emitAcc;
        const nA = Math.max(1, Math.round(nAcross * width / this.d0));
        for (let a = 0; a < nA; a++) {
          const x = this.xN - width / 2 + (a + 0.5) * width / nA + (Math.random() - 0.5) * 0.02 * dp;
          this.addParticle(x, yy, 0, -v, 0, this.yN);
          this.emittedSinceDrain++;
        }
      }
    }

    // 水位保持: 放出した数だけ、コップの底付近の水粒子を取り除く
    drain() {
      if (!this.levelHold || this.emittedSinceDrain <= 0) return;
      const c = this.cup;
      let need = this.emittedSinceDrain;
      const cand = [];
      for (let i = 0; i < this.np && cand.length < need * 6; i++) {
        if (this.state[i] !== 2) continue;
        const [lx, ly] = c.toLocal(this.pos[2 * i], this.pos[2 * i + 1]);
        if (ly > 0 && ly < 0.012 && Math.abs(lx) < c.rb * 0.8) cand.push(i);
      }
      cand.sort((a, b) => b - a); // 後ろから消すと添字がずれない
      let removed = 0;
      for (const i of cand) {
        if (removed >= need) break;
        if (Math.random() < 0.5) continue;
        this.removeParticle(i); removed++;
      }
      this.emittedSinceDrain -= removed;
      if (this.emittedSinceDrain > 2000) this.emittedSinceDrain = 2000;
    }

    // ---------------------------------------------------------- 粒子の積分・衝突
    integrate(dt) {
      const { pos, vel } = this;
      const vcap = 6.0; // 数値的な暴走を防ぐ上限（実際の流れは 3 m/s 以下）
      for (let i = 0; i < this.np; i++) {
        vel[2 * i + 1] -= G * dt;
        const sp = Math.hypot(vel[2 * i], vel[2 * i + 1]);
        if (sp > vcap) { vel[2 * i] *= vcap / sp; vel[2 * i + 1] *= vcap / sp; }
        pos[2 * i] += vel[2 * i] * dt;
        pos[2 * i + 1] += vel[2 * i + 1] * dt;
      }
    }

    pushApart(iters) {
      const { pos, pnx, pny, pInv, cellCount, cellIds } = this;
      const np = this.np;
      cellCount.fill(0);
      for (let i = 0; i < np; i++) {
        const xi = Math.min(Math.max(Math.floor(pos[2 * i] * pInv), 0), pnx - 1);
        const yi = Math.min(Math.max(Math.floor(pos[2 * i + 1] * pInv), 0), pny - 1);
        cellCount[xi * pny + yi]++;
      }
      let first = 0;
      for (let k = 0; k < pnx * pny; k++) { first += cellCount[k]; cellCount[k] = first; }
      cellCount[pnx * pny] = first;
      for (let i = 0; i < np; i++) {
        const xi = Math.min(Math.max(Math.floor(pos[2 * i] * pInv), 0), pnx - 1);
        const yi = Math.min(Math.max(Math.floor(pos[2 * i + 1] * pInv), 0), pny - 1);
        cellIds[--cellCount[xi * pny + yi]] = i;
      }
      const minD = 2 * this.pr, minD2 = minD * minD;
      for (let it = 0; it < iters; it++) {
        for (let i = 0; i < np; i++) {
          const px = pos[2 * i], py = pos[2 * i + 1];
          const pxi = Math.floor(px * pInv), pyi = Math.floor(py * pInv);
          const x0 = Math.max(pxi - 1, 0), y0 = Math.max(pyi - 1, 0);
          const x1 = Math.min(pxi + 1, pnx - 1), y1 = Math.min(pyi + 1, pny - 1);
          for (let xi = x0; xi <= x1; xi++) for (let yi = y0; yi <= y1; yi++) {
            const c = xi * pny + yi;
            for (let q = cellCount[c]; q < cellCount[c + 1]; q++) {
              const id = cellIds[q];
              if (id === i) continue;
              const qx = pos[2 * id], qy = pos[2 * id + 1];
              let dx = qx - px, dy = qy - py;
              const d2 = dx * dx + dy * dy;
              if (d2 > minD2 || d2 === 0) continue;
              const d = Math.sqrt(d2), s = 0.5 * (minD - d) / d;
              dx *= s; dy *= s;
              pos[2 * i] -= dx; pos[2 * i + 1] -= dy;
              pos[2 * id] += dx; pos[2 * id + 1] += dy;
            }
          }
        }
      }
    }

    collide() {
      const { pos, vel, cup, h } = this;
      const minD = cup.rad + this.pr;
      const xmin = h + this.pr, xmax = this.W - h - this.pr, ymin = h + this.pr;
      for (let i = 0; i < this.np; i++) {
        let x = pos[2 * i], y = pos[2 * i + 1];
        // ドメイン外（床・左右）に落ちた水はこぼれた水として消す
        if (x < xmin || x > xmax || y < ymin) {
          if (this.state[i] === 2) this.stats.spilled += this.Vp;
          this.removeParticle(i); i--; continue;
        }
        const bb = cup.bbox;
        const inBox = x > bb[0] && x < bb[1] && y > bb[2] && y < bb[3];
        const dw = inBox ? cup.dist2(x, y) : 1e9;
        // 壁から 2 セル以内を流れる水柱は「壁膜」の一部とみなす（3D では膜に広がる）
        if (this.state[i] === 0 && dw < minD + 2 * h) { this.state[i] = 1; this.yWall[i] = y; this.rough[i] = 0.25 * this.roughAt(this.yb[i] - y); }
        if (dw < minD) {
          const d = cup._d, nx = cup._nx, ny = cup._ny;
          x += nx * (minD - d); y += ny * (minD - d);
          const [wx, wy] = cup.wallVel(x, y);
          const rvx = vel[2 * i] - wx, rvy = vel[2 * i + 1] - wy;
          const vn = rvx * nx + rvy * ny;
          if (vn < 0) { vel[2 * i] -= vn * nx; vel[2 * i + 1] -= vn * ny; }
        }
        // サーバー本体
        if (y > this.yN && Math.abs(x - this.xN) < 0.03 && !(Math.abs(x - this.xN) < this.d0 / 2 && this.state[i] === 0)) {
          if (vel[2 * i + 1] > 0) {
            // 上に飛んだ大粒がノズルに当たった
            if (Math.abs(x - this.xN) < this.rNoz * 1.6 && this.state[i] === 2) this.stats.bigNozzleHits += 1;
            vel[2 * i + 1] = 0;
          }
          y = this.yN;
        }
        pos[2 * i] = x; pos[2 * i + 1] = y;
      }
    }

    roughAt(zfall) {
      return Math.min(this.noise * Math.exp(interpTable(this.tab, this.tab.gain, Math.max(zfall, 0))), 1);
    }

    // ---------------------------------------------------------- 粒子 ↔ 格子
    transfer(toGrid) {
      const { nx, ny, h, pos, vel, u, v, du, dv, pu, pv, s, type } = this;
      const h1 = 1 / h, h2 = 0.5 * h;
      const np = this.np;
      if (toGrid) {
        u.fill(0); v.fill(0); du.fill(0); dv.fill(0);
        for (let k = 0; k < nx * ny; k++) type[k] = s[k] === 0 ? SOLID : AIR;
        for (let i = 0; i < np; i++) {
          const xi = Math.min(Math.max(Math.floor(pos[2 * i] * h1), 0), nx - 1);
          const yi = Math.min(Math.max(Math.floor(pos[2 * i + 1] * h1), 0), ny - 1);
          const k = xi * ny + yi;
          if (type[k] === AIR) type[k] = FLUID;
        }
      }
      for (let comp = 0; comp < 2; comp++) {
        const dx = comp === 0 ? 0 : h2, dy = comp === 0 ? h2 : 0;
        const f = comp === 0 ? u : v, prevF = comp === 0 ? pu : pv, d = comp === 0 ? du : dv;
        for (let i = 0; i < np; i++) {
          let x = pos[2 * i], y = pos[2 * i + 1];
          x = Math.min(Math.max(x, h), (nx - 1) * h);
          y = Math.min(Math.max(y, h), (ny - 1) * h);
          const x0 = Math.min(Math.floor((x - dx) * h1), nx - 2), tx = ((x - dx) - x0 * h) * h1, x1 = Math.min(x0 + 1, nx - 2);
          const y0 = Math.min(Math.floor((y - dy) * h1), ny - 2), ty = ((y - dy) - y0 * h) * h1, y1 = Math.min(y0 + 1, ny - 2);
          const sx = 1 - tx, sy = 1 - ty;
          const d0 = sx * sy, d1 = tx * sy, d2 = tx * ty, d3 = sx * ty;
          const n0 = x0 * ny + y0, n1 = x1 * ny + y0, n2 = x1 * ny + y1, n3 = x0 * ny + y1;
          if (toGrid) {
            const pvv = vel[2 * i + comp];
            f[n0] += pvv * d0; d[n0] += d0;
            f[n1] += pvv * d1; d[n1] += d1;
            f[n2] += pvv * d2; d[n2] += d2;
            f[n3] += pvv * d3; d[n3] += d3;
          } else {
            // 固体に接する面（階段状の壁で接線成分まで 0 になる）は使わない → 斜めの壁でも膜が滑る
            const off = comp === 0 ? ny : 1;
            const ok = (a, b) => (a !== SOLID && b !== SOLID && (a === FLUID || b === FLUID)) ? 1 : 0;
            const v0 = ok(type[n0], type[n0 - off]);
            const v1 = ok(type[n1], type[n1 - off]);
            const v2 = ok(type[n2], type[n2 - off]);
            const v3 = ok(type[n3], type[n3 - off]);
            const sw = v0 * d0 + v1 * d1 + v2 * d2 + v3 * d3;
            if (sw > 0) {
              const pic = (v0 * d0 * f[n0] + v1 * d1 * f[n1] + v2 * d2 * f[n2] + v3 * d3 * f[n3]) / sw;
              const corr = (v0 * d0 * (f[n0] - prevF[n0]) + v1 * d1 * (f[n1] - prevF[n1])
                + v2 * d2 * (f[n2] - prevF[n2]) + v3 * d3 * (f[n3] - prevF[n3])) / sw;
              const flipV = vel[2 * i + comp] + corr;
              vel[2 * i + comp] = (1 - this.flip) * pic + this.flip * flipV;
            }
          }
        }
        if (toGrid) {
          for (let k = 0; k < f.length; k++) if (d[k] > 0) f[k] /= d[k];
          // 固体に接する面は固体の速度
          for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
            const k = i * ny + j;
            const solid = type[k] === SOLID;
            if (comp === 0 && i > 0 && (solid || type[k - ny] === SOLID)) {
              const kk = solid ? k : k - ny;
              u[k] = this.solidVx[kk];
            }
            if (comp === 1 && j > 0 && (solid || type[k - 1] === SOLID)) {
              const kk = solid ? k : k - 1;
              v[k] = this.solidVy[kk];
            }
          }
        }
      }
    }

    updateDensity() {
      const { nx, ny, h, pos, dens, poolCnt } = this;
      const h1 = 1 / h, h2 = 0.5 * h;
      dens.fill(0); poolCnt.fill(0); this.filmCnt.fill(0);
      for (let i = 0; i < this.np; i++) {
        let x = Math.min(Math.max(pos[2 * i], h), (nx - 1) * h);
        let y = Math.min(Math.max(pos[2 * i + 1], h), (ny - 1) * h);
        const x0 = Math.floor((x - h2) * h1), tx = ((x - h2) - x0 * h) * h1, x1 = Math.min(x0 + 1, nx - 2);
        const y0 = Math.floor((y - h2) * h1), ty = ((y - h2) - y0 * h) * h1, y1 = Math.min(y0 + 1, ny - 2);
        const sx = 1 - tx, sy = 1 - ty;
        dens[x0 * ny + y0] += sx * sy; dens[x1 * ny + y0] += tx * sy;
        dens[x1 * ny + y1] += tx * ty; dens[x0 * ny + y1] += sx * ty;
        if (this.state[i] === 2) {
          const xi = Math.floor(x * h1), yi = Math.floor(y * h1);
          poolCnt[xi * ny + yi] += 1;
        } else if (this.state[i] === 1) {
          const xi = Math.floor(x * h1), yi = Math.floor(y * h1);
          this.filmCnt[xi * ny + yi] += 1;
        }
      }
    }

    solvePressure(dt) {
      const { nx, ny, u, v, s, type, dens } = this;
      const omega = 1.9;
      const rest = this.restDensity;
      // FLIP の差分は「圧力で変わった分」だけにする（P2G 直後の速度を保存）
      this.pu.set(u); this.pv.set(v);
      // 水セルの一覧を作ってから反復（空気セルを毎回なめない）
      if (!this._fl || this._fl.length < nx * ny) this._fl = new Int32Array(nx * ny);
      const fl = this._fl;
      let nf = 0;
      for (let i = 1; i < nx - 1; i++) for (let j = 1; j < ny - 1; j++) {
        const k = i * ny + j;
        if (type[k] === FLUID) fl[nf++] = k;
      }
      for (let it = 0; it < this.pressureIters; it++) {
        for (let q = 0; q < nf; q++) {
          const k = fl[q];
          const sx0 = s[k - ny], sx1 = s[k + ny], sy0 = s[k - 1], sy1 = s[k + 1];
          const st = sx0 + sx1 + sy0 + sy1;
          if (st === 0) continue;
          let div = u[k + ny] - u[k] + v[k + 1] - v[k];
          if (rest > 0) {
            const comp = dens[k] - rest;
            if (comp > 0) div -= 1.0 * comp;
          }
          const pp = -div / st * omega;
          u[k] -= sx0 * pp; u[k + ny] += sx1 * pp;
          v[k] -= sy0 * pp; v[k + 1] += sy1 * pp;
        }
      }
    }

    // ---------------------------------------------------------- サブグリッド: 着水判定
    detectImpacts() {
      const { pos, vel, h, ny } = this;
      const h1 = 1 / h;
      const st = this.stats;
      for (let i = 0; i < this.np; i++) {
        let s0 = this.state[i];
        if (s0 === 2) continue;
        const x = pos[2 * i], y = pos[2 * i + 1];
        const k = Math.floor(x * h1) * ny + Math.floor(y * h1);
        // 壁膜に合流した水柱の粒子も壁膜扱い（膜の外側の層）
        if (s0 === 0) {
          const fc = this.filmCnt;
          // 条件: 膜の粒子が隣にいて、壁から 5 セル以内で、速度が壁にほぼ平行（壁で曲げられた後）
          if (fc[k] + fc[k - 1] + fc[k + 1] + fc[k - ny] + fc[k + ny] > 0
              && this.cup.dist2(x, y) < this.cup.rad + this.pr + 5 * h
              && Math.abs(vel[2 * i] * this.cup._nx + vel[2 * i + 1] * this.cup._ny) < 0.35 * Math.hypot(vel[2 * i], vel[2 * i + 1])) {
            this.state[i] = 1; s0 = 1; this.yWall[i] = y; this.rough[i] = 0.25 * this.roughAt(this.yb[i] - y);
          }
        }
        const below = k - 1;
        // 飛び散った粒ではなく「水のかたまり」に入ったときだけ着水とみなす
        if (this.poolCnt[k] + 0.5 * this.poolCnt[below] < 2.0 || this.dens[k] < 0.5 * this.restDensity) continue;
        // 着水した。速さは格子の数値粘性で鈍るので、落下距離から理論値で評価する（3D の水柱は空中で減速しない）
        const vex = Math.max(this.exitVelocity(), 0.3 * this.v0);
        let Hf, U, Ueff, rough;
        if (s0 === 0) {
          Hf = Math.max(this.yN - y, 1e-3);
          U = Math.sqrt(vex * vex + 2 * G * Hf);
          Ueff = U;
          rough = this.roughAt(Hf);
        } else {
          // 壁膜: 壁に当たった高さ・速さ → 垂直成分の半分を失い → 壁沿いに重力で加速 → 水面に垂直な成分だけ突っ込む
          Hf = Math.max(this.yN - this.yWall[i], 1e-3);
          const Uw = Math.sqrt(vex * vex + 2 * G * Hf);
          const phi = Math.abs(this.cup.th) + Math.atan((this.cup.rt - this.cup.rb) / this.cup.hc);
          const Un = Uw * Math.sin(phi);
          const drop = Math.max(this.yWall[i] - y, 0);
          const uf = Math.sqrt(Math.max(Uw * Uw - 0.5 * Un * Un + 2 * G * drop, 0));
          U = Uw;
          Ueff = uf * Math.cos(phi);
          rough = this.rough[i];
        }
        const dJet = 2 * interpTable(this.tab, this.tab.a, Hf);
        let ar = airRatio(Ueff, Hf, this.d0, this.v0, rough, this.fl);
        if (s0 === 1) ar *= 0.3;
        const Ve = onsetVelocity(rough, this.fl);
        // 指数移動平均で表示値を更新
        const a = 0.02;
        st.H += a * (Hf - st.H);
        st.U += a * (U - st.U);
        st.Ueff += a * (Ueff - st.Ueff);
        st.rough += a * (rough - st.rough);
        st.Ve += a * (Ve - st.Ve);
        st.air += a * (ar - st.air);
        st.wallFrac += a * ((s0 === 1 ? 1 : 0) - st.wallFrac);
        st.hasImpact = true;
        this.state[i] = 2;
        // 巻き込んだ空気 → 気泡のかたまり（1 個のパーセルが同じ大きさの n 個の気泡を代表）
        const Va = ar * this.Vp;
        if (Va > 0) {
          const Rmed = bubbleRmed(Ueff, Math.max(dJet, 1e-4), this.fl);
          const R = Math.min(Math.max(Rmed * Math.exp(0.6 * randn()), 20e-6), 2.5e-3);
          const nreal = Va / (4 / 3 * Math.PI * R * R * R);
          this.spawnBubble(x + (Math.random() - 0.5) * dJet, y - Math.random() * 0.004, R, nreal);
        }
      }
    }

    spawnBubble(x, y, R, w) {
      this.stats.bubbles += w; this._accB += w;
      if (this.nb < this.maxB) {
        const k = this.nb++;
        this.bx[k] = x; this.by[k] = y; this.bR[k] = R; this.bw[k] = w;
      } else {
        const k = Math.floor(Math.random() * this.nb);
        this.bw[k] += w; // 表示上の上限を超えたら既存のパーセルに重みを足す
      }
    }

    sampleGrid(x, y, comp) {
      const { nx, ny, h } = this;
      const h1 = 1 / h, h2 = 0.5 * h, f = comp === 0 ? this.u : this.v;
      const dx = comp === 0 ? 0 : h2, dy = comp === 0 ? h2 : 0;
      x = Math.min(Math.max(x, h), (nx - 1) * h); y = Math.min(Math.max(y, h), (ny - 1) * h);
      const x0 = Math.min(Math.floor((x - dx) * h1), nx - 2), tx = ((x - dx) - x0 * h) * h1;
      const y0 = Math.min(Math.floor((y - dy) * h1), ny - 2), ty = ((y - dy) - y0 * h) * h1;
      return (1 - tx) * (1 - ty) * f[x0 * ny + y0] + tx * (1 - ty) * f[(x0 + 1) * ny + y0]
        + tx * ty * f[(x0 + 1) * ny + y0 + 1] + (1 - tx) * ty * f[x0 * ny + y0 + 1];
    }
    cellDensity(x, y) {
      const xi = Math.floor(x / this.h), yi = Math.floor(y / this.h);
      if (xi < 0 || yi < 0 || xi >= this.nx || yi >= this.ny) return 0;
      return this.dens[xi * this.ny + yi] / Math.max(this.restDensity, 1e-6);
    }
    cellSolid(x, y) {
      const xi = Math.floor(x / this.h), yi = Math.floor(y / this.h);
      if (xi < 0 || yi < 0 || xi >= this.nx || yi >= this.ny) return true;
      return this.s[xi * this.ny + yi] === 0;
    }

    stepBubbles(dt) {
      for (let k = 0; k < this.nb; k++) {
        const x = this.bx[k], y = this.by[k], R = this.bR[k];
        const ux = this.sampleGrid(x, y, 0), uy = this.sampleGrid(x, y, 1);
        const vr = riseVelocity(R, this.fl);
        let nx = x + ux * dt, ny = y + (uy + vr) * dt;
        if (this.cellSolid(nx, ny)) { nx = x; ny = y + vr * dt; }
        this.bx[k] = nx; this.by[k] = ny;
        const rho = this.cellDensity(nx, ny);
        if (rho < 0.45 || this.cellDensity(nx, ny + this.h) < 0.12) {
          this.burst(k);
          // 末尾と入れ替えて削除
          const last = --this.nb;
          this.bx[k] = this.bx[last]; this.by[k] = this.by[last]; this.bR[k] = this.bR[last]; this.bw[k] = this.bw[last];
          k--;
        }
      }
    }

    burst(k) {
      const R = this.bR[k], w = this.bw[k];
      const jd = jetDrops(R, this.fl);
      for (let m = 0; m < jd.n; m++) {
        const r = jd.rd * Math.pow(1.08, m), V = jd.Vd * Math.pow(0.6, m);
        const ang = Math.PI / 2 + 0.14 * randn();
        this.spawnDrop(this.bx[k], this.by[k] + 0.0005, V * Math.cos(ang), V * Math.sin(ang), r, w);
      }
    }

    spawnDrop(x, y, vx, vy, r, w) {
      this.stats.drops += w; this._accD += w;
      let k;
      if (this.nd < this.maxD) k = this.nd++;
      else {
        k = Math.floor(Math.random() * this.nd);
        // 押し出される古いパーセルの重みは捨てる（統計は spawn 時に加算済み）
      }
      this.dx_[k] = x; this.dy_[k] = y; this.dvx[k] = vx; this.dvy[k] = vy; this.dr[k] = r; this.dw[k] = w;
      this.dTop[k] = y; this.dFlag[k] = 0;
    }

    stepDrops(dt) {
      const rhoA = 1.184, muA = 1.85e-5, rhoL = this.fl.rho;
      const rimTop = Math.max(this.cup.rimL[1], this.cup.rimR[1]);
      for (let k = 0; k < this.nd; k++) {
        const r = this.dr[k];
        let vx = this.dvx[k], vy = this.dvy[k];
        const vm = Math.hypot(vx, vy) + 1e-9;
        const Re = 2 * r * vm * rhoA / muA;
        const Cd = Re < 1000 ? 24 / Math.max(Re, 1e-9) * (1 + 0.15 * Math.pow(Re, 0.687)) : 0.44;
        const kd = 3 * rhoA * Cd * vm / (8 * rhoL * r);   // 抵抗の緩和率 [1/s]
        // 指数積分（小さい粒でも安定）: dv/dt = -kd v + g
        const e = Math.exp(-kd * dt), vt = -G / kd;
        vx = vx * e;
        vy = vt + (vy - vt) * e;
        const x = this.dx_[k] + vx * dt, y = this.dy_[k] + vy * dt;
        this.dvx[k] = vx; this.dvy[k] = vy; this.dx_[k] = x; this.dy_[k] = y;
        if (y > this.dTop[k]) this.dTop[k] = y;
        let dead = false;
        // ノズルに当たった
        if (y >= this.yN - 0.0005 && Math.abs(x - this.xN) < this.rNoz) {
          this.stats.nozzleHits += this.dw[k]; dead = true;
          this.stats.events.push({ t: this.time, x, y });
        } else if (y > rimTop && !(this.dFlag[k] & 1)) {
          this.dFlag[k] |= 1;
          this.stats.rimEscape += this.dw[k];
        }
        if (!dead && vy < 0 && (this.cellDensity(x, y) > 0.35 || this.cellSolid(x, y))) dead = true;
        if (!dead && (x < 0 || x > this.W || y < 0 || y > this.H)) dead = true;
        if (dead) {
          const last = --this.nd;
          this.dx_[k] = this.dx_[last]; this.dy_[k] = this.dy_[last]; this.dvx[k] = this.dvx[last]; this.dvy[k] = this.dvy[last];
          this.dr[k] = this.dr[last]; this.dw[k] = this.dw[last]; this.dTop[k] = this.dTop[last]; this.dFlag[k] = this.dFlag[last];
          k--;
        }
      }
    }

    // ---------------------------------------------------------- 1 ステップ
    substep(dt) {
      this.emit(dt);
      this.drain();
      this.integrate(dt);
      this.pushApart(1);
      this.collide();
      this.transfer(true);
      this.updateDensity();
      this.solvePressure(dt);
      this.transfer(false);
      this.detectImpacts();
      this.stepBubbles(dt);
      this.stepDrops(dt);
      this.time += dt;
      // 発生率（1 秒窓）
      this._accT += dt;
      if (this._accT > 0.25) {
        const st = this.stats;
        st.bubbleRate = this._accB / this._accT; st.dropRate = this._accD / this._accT;
        this._accB = 0; this._accD = 0; this._accT = 0;
      }
    }

    step(frameDt) {
      // CFL で部分ステップ数を決める
      let vmax = 0.5;
      for (let i = 0; i < this.np; i++) {
        const a = Math.abs(this.vel[2 * i]), b = Math.abs(this.vel[2 * i + 1]);
        if (a > vmax) vmax = a; if (b > vmax) vmax = b;
      }
      vmax = Math.max(vmax, 1.0);
      const dtMax = 1.2 * this.h / vmax;
      const n = Math.min(Math.max(1, Math.ceil(frameDt / dtMax)), 40);
      const dt = frameDt / n;
      for (let k = 0; k < n; k++) this.substep(dt);
      return n;
    }

    // コップを中の水ごと瞬間移動させる（初期配置・プリセット用）
    teleportCup(x, y, th) {
      const c = this.cup;
      const ids = [], loc = [];
      for (let i = 0; i < this.np; i++) {
        const px = this.pos[2 * i], py = this.pos[2 * i + 1];
        const [lx, ly] = c.toLocal(px, py);
        const hw = c.rb + (c.rt - c.rb) * Math.min(Math.max(ly, 0), c.hc) / c.hc + 0.003;
        if (ly > -0.003 && ly < c.hc + 0.01 && Math.abs(lx) < hw) { ids.push(i); loc.push([lx, ly]); }
      }
      c.x = x; c.y = y; c.th = th; c.vx = c.vy = c.om = 0;
      c.update();
      ids.forEach((i, k) => {
        const [wx, wy] = c.toWorld(loc[k][0], loc[k][1]);
        this.pos[2 * i] = wx; this.pos[2 * i + 1] = wy;
        this.vel[2 * i] = 0; this.vel[2 * i + 1] = 0;
      });
      // 傾けて水面より上に出た分（こぼれる分）はそのまま流れ落ちる
      this.u.fill(0); this.v.fill(0);
      this.markSolids();
    }

    // コップの姿勢を変える（壁の速度も計算）
    moveCup(x, y, th, dt) {
      const c = this.cup;
      if (dt > 0) {
        c.vx = (x - c.x) / dt; c.vy = (y - c.y) / dt;
        let dth = th - c.th;
        c.om = dth / dt;
        const vmax = 1.2;
        const sp = Math.hypot(c.vx, c.vy);
        if (sp > vmax) { c.vx *= vmax / sp; c.vy *= vmax / sp; }
      } else { c.vx = c.vy = c.om = 0; }
      c.x = x; c.y = y; c.th = th;
      this.markSolids();
    }

    // ---------------------------------------------------------- 注ぎ方 → コップの目標姿勢（解析的）
    // 傾き th のコップに 2D 面積 A の水が入っているときの、内底中心から測った水面の高さ（世界 y の差）
    surfaceOffset(A, th) {
      const c = this.cup, cs = Math.cos(th), sn = Math.sin(th);
      const loc = [[-c.rb, 0], [c.rb, 0], [c.rt, c.hc], [-c.rt, c.hc]];
      const poly = loc.map(([x, y]) => [x * cs - y * sn, x * sn + y * cs]);
      const areaBelow = (yc) => {
        const out = [];
        for (let i = 0; i < 4; i++) {
          const a = poly[i], b = poly[(i + 1) % 4];
          const ain = a[1] <= yc, bin = b[1] <= yc;
          if (ain) out.push(a);
          if (ain !== bin) { const t = (yc - a[1]) / (b[1] - a[1]); out.push([a[0] + t * (b[0] - a[0]), yc]); }
        }
        let ar = 0;
        for (let i = 0; i < out.length; i++) { const a = out[i], b = out[(i + 1) % out.length]; ar += a[0] * b[1] - b[0] * a[1]; }
        return Math.abs(ar) / 2;
      };
      let lo = Math.min(...poly.map(p => p[1])), hi = Math.max(...poly.map(p => p[1]));
      for (let it = 0; it < 50; it++) { const m = 0.5 * (lo + hi); if (areaBelow(m) < A) lo = m; else hi = m; }
      const rimLow = Math.min(poly[2][1], poly[3][1]);
      return { ys: 0.5 * (lo + hi), rimLow };
    }
    waterArea() {
      let n = 0;
      for (let i = 0; i < this.np; i++) if (this.state[i] === 2 && this.cup.inside(this.pos[2 * i], this.pos[2 * i + 1])) n++;
      return n * this.dp * this.dp;
    }
    // 戦略からコップの目標姿勢を出す。aim='wall' で着地点が縁に近すぎるときは傾きを自動で浅くする
    targetPose(gap, thDeg, aim = 'center', offset = 0.02, A = null) {
      const c = this.cup, k = (c.rt - c.rb) / c.hc;
      A = A ?? this.waterArea();
      const r = c.rt + c.tw / 2;
      let deg = thDeg;
      for (let tries = 0; tries < 40; tries++) {
        const th = deg * Math.PI / 180, cs = Math.cos(th), sn = Math.sin(th);
        const { ys, rimLow } = this.surfaceOffset(A, th);
        const rimTop = Math.max(-r * sn + c.hc * cs, r * sn + c.hc * cs);
        let px, ok = ys < rimLow - 0.006;
        if (aim === 'wall') {
          const side = th >= 0 ? -1 : 1;
          // 床側内壁 lx = side (rb + k ly) の、世界 y が ys + offset になる点
          let ly = (ys + offset - side * c.rb * sn) / (side * k * sn + cs);
          if (ly > c.hc - 0.018) ok = false;
          ly = Math.min(Math.max(ly, 0.004), c.hc - 0.012);
          const lx = side * (c.rb + k * ly) + (Math.abs(th) < 1e-6 ? 0.0012 : 0);
          px = lx * cs - ly * sn;
        } else {
          const xs = [];
          for (const sd of [-1, 1]) {
            const l = (ys - sd * c.rb * sn) / (sd * k * sn + cs);
            xs.push(sd * (c.rb + k * l) * cs - l * sn);
          }
          px = 0.5 * (xs[0] + xs[1]);
        }
        if (ok || Math.abs(deg) < 1) return { x: this.xN - px, y: this.yN - gap - rimTop, th, deg, ok };
        deg = deg > 0 ? Math.max(deg - 2, 0) : Math.min(deg + 2, 0);
      }
      return null;
    }

    // コップの中の水の割合（2D の面積比）
    cupFill() {
      let n = 0;
      for (let i = 0; i < this.np; i++) if (this.state[i] === 2 && this.cup.inside(this.pos[2 * i], this.pos[2 * i + 1])) n++;
      return n * this.dp * this.dp / this.cup.capacityArea();
    }
    // ノズル直下の水面の高さ（ノズルから水面までの距離）
    surfaceBelowNozzle() {
      const i = Math.floor(this.xN / this.h);
      const thr = 0.5 * this.restDensity;
      for (let j = Math.floor(this.yN / this.h) - 2; j > 0; j--) {
        const k = i * this.ny + j;
        if (this.s[k] === 0) return { y: (j + 1) * this.h, kind: 'wall' };
        if (this.poolCnt[k] >= 2 && this.dens[k] > thr) return { y: (j + 1) * this.h, kind: 'water' };
      }
      return null;
    }
  }

  const api = { PourSim, Cup, waterProps, jetTable, onsetVelocity, airRatio, bubbleRmed, jetDrops, G };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.PourPhys = api;
})(typeof window !== 'undefined' ? window : globalThis);
