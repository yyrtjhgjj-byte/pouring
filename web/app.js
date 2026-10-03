/*
 * app.js — 注水シミュレーターの画面（描画と操作）
 * sim.js の PourSim を動かし、キャンバスに描いて、計測値を表示する。
 */
(function () {
  'use strict';
  const { PourSim } = window.PourPhys;
  const $ = (id) => document.getElementById(id);
  const canvas = $('sim-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  // ------------------------------------------------------------------ 状態
  const ui = {
    gap: 0.12, tilt: 0, aim: 'center', autoAim: true,
    speed: 0.25, latched: true, holding: false, showStates: false,
    playing: true,
  };
  let sim;
  let pose = null;            // 現在のコップ姿勢 {x,y,th}
  let flashUntil = 0, lastHits = 0;
  let dragging = null;

  const PRESETS = {
    far: { gap: 0.12, tilt: 0, aim: 'center', tau: 0.06, label: '遠くから中央へ' },
    near: { gap: 0.01, tilt: 0, aim: 'center', tau: 0.06, label: '近づけて中央へ' },
    best: { gap: 0.025, tilt: 30, aim: 'wall', tau: 0.25, label: '傾けて壁沿い' },
  };

  function newSim() {
    const T = +document.querySelector('input[name="temp"]:checked').value;
    sim = new PourSim({ h: 0.00125, fillFrac: 0.3, T, Q: +$('flow').value / 60000, tauValve: +$('tau').value, pressureIters: 24 });
    sim.levelHold = $('level-hold').checked;
    const tp = sim.targetPose(ui.gap, ui.tilt, ui.aim, 0.012);
    sim.teleportCup(tp.x, tp.y, tp.th);
    pose = { x: tp.x, y: tp.y, th: tp.th };
    sim.lever = false;
    lastHits = 0;
  }

  // ------------------------------------------------------------------ 操作
  function syncLabels() {
    $('gap-val').textContent = (ui.gap * 100).toFixed(1) + ' cm';
    $('tilt-val').textContent = ui.tilt.toFixed(0) + '°';
    $('tau-val').textContent = (+$('tau').value).toFixed(2) + ' s';
    $('flow-val').textContent = (+$('flow').value).toFixed(1) + ' L/分';
    $('gap').value = ui.gap * 100;
    $('tilt').value = ui.tilt;
    document.querySelectorAll('input[name="aim"]').forEach(r => { r.checked = r.value === ui.aim; });
    $('auto-aim').checked = ui.autoAim;
  }

  function applyPreset(key) {
    const p = PRESETS[key];
    ui.gap = p.gap; ui.tilt = p.tilt; ui.aim = p.aim; ui.autoAim = true;
    $('tau').value = p.tau; sim.tauValve = p.tau;
    document.querySelectorAll('.preset').forEach(b => b.setAttribute('aria-pressed', b.dataset.preset === key ? 'true' : 'false'));
    // いったん止めて、コップを水ごと目標の姿勢に置き、水面が落ち着いてから注ぎ直す（開栓の瞬間も見えるように）
    sim.lever = false;
    const tp = sim.targetPose(ui.gap, ui.tilt, ui.aim, 0.012);
    if (tp) { sim.teleportCup(tp.x, tp.y, tp.th); pose = { x: tp.x, y: tp.y, th: tp.th }; }
    pendingPour = 0.5; // シミュレーション時間でこの後に注ぎ始める
    syncLabels(); updatePourBtn();
  }
  let pendingPour = 0.4;

  $('gap').addEventListener('input', e => { ui.gap = e.target.value / 100; syncLabels(); clearPreset(); });
  $('tilt').addEventListener('input', e => { ui.tilt = +e.target.value; syncLabels(); clearPreset(); });
  document.querySelectorAll('input[name="aim"]').forEach(r => r.addEventListener('change', e => { ui.aim = e.target.value; ui.autoAim = true; syncLabels(); clearPreset(); }));
  $('auto-aim').addEventListener('change', e => { ui.autoAim = e.target.checked; });
  $('tau').addEventListener('input', e => { sim.tauValve = +e.target.value; syncLabels(); });
  $('flow').addEventListener('input', e => { sim.setFlow(+e.target.value / 60000); syncLabels(); });
  document.querySelectorAll('input[name="temp"]').forEach(r => r.addEventListener('change', e => sim.setTemperature(+e.target.value)));
  document.querySelectorAll('input[name="speed"]').forEach(r => r.addEventListener('change', e => { ui.speed = +e.target.value; }));
  $('level-hold').addEventListener('change', e => { sim.levelHold = e.target.checked; });
  $('show-states').addEventListener('change', e => { ui.showStates = e.target.checked; });
  $('latch').addEventListener('change', e => { ui.latched = e.target.checked; if (!ui.latched) sim.lever = false; updatePourBtn(); });
  document.querySelectorAll('.preset').forEach(b => b.addEventListener('click', () => applyPreset(b.dataset.preset)));
  $('reset').addEventListener('click', () => { newSim(); pendingPour = 0.4; });
  $('clear-counts').addEventListener('click', () => { sim.resetStats(); lastHits = 0; });
  $('play').addEventListener('click', () => { ui.playing = !ui.playing; $('play').textContent = ui.playing ? '一時停止' : '再生'; if (ui.playing) last = performance.now(); });
  function clearPreset() { document.querySelectorAll('.preset').forEach(b => b.setAttribute('aria-pressed', 'false')); }

  const pourBtn = $('pour');
  function updatePourBtn() {
    pourBtn.setAttribute('aria-pressed', sim && sim.lever ? 'true' : 'false');
    pourBtn.textContent = sim && sim.lever ? (ui.latched ? '止める' : '注いでいる…') : (ui.latched ? '注ぐ' : '押している間だけ注ぐ');
  }
  pourBtn.addEventListener('pointerdown', (e) => {
    e.preventDefault();
    if (ui.latched) { sim.lever = !sim.lever; pendingPour = -1; }
    else { sim.lever = true; ui.holding = true; pendingPour = -1; }
    updatePourBtn();
  });
  const release = () => { if (!ui.latched && ui.holding) { ui.holding = false; sim.lever = false; updatePourBtn(); } };
  pourBtn.addEventListener('pointerup', release);
  pourBtn.addEventListener('pointerleave', release);
  window.addEventListener('keydown', (e) => {
    if (e.code === 'Space' && e.target === document.body) { e.preventDefault(); if (ui.latched) sim.lever = !sim.lever; else sim.lever = true; pendingPour = -1; updatePourBtn(); }
  });
  window.addEventListener('keyup', (e) => { if (e.code === 'Space' && !ui.latched) { sim.lever = false; updatePourBtn(); } });

  // コップをドラッグ
  function toWorld(ev) {
    const r = canvas.getBoundingClientRect();
    const sx = (ev.clientX - r.left) / r.width, sy = (ev.clientY - r.top) / r.height;
    return [sx * sim.W, (1 - sy) * sim.H];
  }
  canvas.addEventListener('pointerdown', (ev) => {
    const [x, y] = toWorld(ev);
    const [lx, ly] = sim.cup.toLocal(x, y);
    if (ly > -0.01 && ly < sim.cup.hc + 0.015 && Math.abs(lx) < sim.cup.rt + 0.012) {
      dragging = { dx: pose.x - x, dy: pose.y - y };
      ui.autoAim = false; syncLabels(); clearPreset();
      canvas.setPointerCapture(ev.pointerId);
    }
  });
  canvas.addEventListener('pointermove', (ev) => {
    if (!dragging) return;
    const [x, y] = toWorld(ev);
    dragTarget = { x: x + dragging.dx, y: y + dragging.dy };
  });
  canvas.addEventListener('pointerup', () => { dragging = null; });
  let dragTarget = null;

  // ------------------------------------------------------------------ コップの動き（目標へなめらかに）
  function moveCupToward(dtSim) {
    let target;
    if (dragTarget && !ui.autoAim) {
      target = { x: dragTarget.x, y: dragTarget.y, th: ui.tilt * Math.PI / 180 };
    } else if (ui.autoAim) {
      const tp = sim.targetPose(ui.gap, ui.tilt, ui.aim, 0.012);
      target = tp ? { x: tp.x, y: tp.y, th: tp.th } : pose;
      if (tp && Math.abs(tp.deg - ui.tilt) > 0.5) { $('tilt-note').textContent = `水が多いので ${tp.deg.toFixed(0)}° に抑えています`; }
      else $('tilt-note').textContent = '';
    } else {
      target = { x: pose.x, y: pose.y, th: ui.tilt * Math.PI / 180 };
    }
    // 床やサーバーにめり込まないように
    target.y = Math.max(target.y, 0.012);
    const vmax = 0.12, wmax = 0.5; // m/s, rad/s（シミュレーション時間）。速く動かすと水がこぼれる
    const dx = target.x - pose.x, dy = target.y - pose.y, dth = target.th - pose.th;
    const d = Math.hypot(dx, dy), stepL = vmax * dtSim, stepA = wmax * dtSim;
    const f = d > stepL ? stepL / d : 1;
    const nth = pose.th + Math.max(-stepA, Math.min(stepA, dth));
    const nx = pose.x + dx * f, ny = pose.y + dy * f;
    const moved = Math.abs(nx - pose.x) + Math.abs(ny - pose.y) + Math.abs(nth - pose.th) > 1e-7;
    pose = { x: nx, y: ny, th: nth };
    if (moved) sim.moveCup(nx, ny, nth, dtSim);
    else if (sim.cup.vx || sim.cup.vy || sim.cup.om) sim.moveCup(nx, ny, nth, 0);
  }

  // ------------------------------------------------------------------ 描画
  let off = document.createElement('canvas'), offCtx = off.getContext('2d'), img = null, field = null;
  function resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(w * dpr * sim.H / sim.W);
  }

  function renderWater(scale) {
    const R = 2; // 格子の 2 倍の細かさで粒子を塗る
    const nx = sim.nx * R, ny = sim.ny * R, hh = sim.h / R;
    if (!field || field.length !== nx * ny) {
      field = new Float32Array(nx * ny);
      off.width = nx; off.height = ny;
      img = offCtx.createImageData(nx, ny);
    }
    field.fill(0);
    const pos = sim.pos, inv = 1 / hh;
    for (let i = 0; i < sim.np; i++) {
      const x = pos[2 * i] * inv - 0.5, y = pos[2 * i + 1] * inv - 0.5;
      const x0 = Math.floor(x), y0 = Math.floor(y);
      if (x0 < 1 || y0 < 1 || x0 >= nx - 2 || y0 >= ny - 2) continue;
      const tx = x - x0, ty = y - y0;
      // 3x3 のやわらかい重み
      for (let a = -1; a <= 2; a++) {
        const wx = a === -1 ? (1 - tx) * 0.25 : a === 2 ? tx * 0.25 : (a === 0 ? 1 - tx * 0.5 : 0.5 + tx * 0.5) * 0.75;
        for (let b = -1; b <= 2; b++) {
          const wy = b === -1 ? (1 - ty) * 0.25 : b === 2 ? ty * 0.25 : (b === 0 ? 1 - ty * 0.5 : 0.5 + ty * 0.5) * 0.75;
          field[(x0 + a) * ny + (y0 + b)] += wx * wy;
        }
      }
    }
    const data = img.data;
    const showS = ui.showStates;
    for (let i = 0; i < nx; i++) {
      for (let j = 0; j < ny; j++) {
        const v = field[i * ny + j];
        const p = ((ny - 1 - j) * nx + i) * 4;
        if (v < 0.12) { data[p + 3] = 0; continue; }
        const a = Math.min(1, (v - 0.12) / 0.28);
        // 表面ほど明るく
        const edge = Math.max(0, 1 - (v - 0.12) / 0.9);
        data[p] = 40 + 110 * edge; data[p + 1] = 130 + 90 * edge; data[p + 2] = 210 + 40 * edge;
        data[p + 3] = 235 * a;
      }
    }
    offCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
    if (showS) {
      // 粒子の状態を色分け: 緑=空中の水柱, 橙=壁を伝う膜
      for (let i = 0; i < sim.np; i++) {
        const st = sim.state[i];
        if (st === 2) continue;
        ctx.fillStyle = st === 0 ? 'rgba(120,230,140,0.9)' : 'rgba(255,150,70,0.95)';
        ctx.fillRect(sim.pos[2 * i] * scale - 1, (sim.H - sim.pos[2 * i + 1]) * scale - 1, 2, 2);
      }
    }
  }

  function wy(y, scale) { return (sim.H - y) * scale; }

  function draw(now) {
    const W = canvas.width, H = canvas.height, scale = W / sim.W;
    // 背景
    const g = ctx.createLinearGradient(0, 0, 0, H);
    g.addColorStop(0, '#0d1b2a'); g.addColorStop(1, '#13283d');
    ctx.fillStyle = g; ctx.fillRect(0, 0, W, H);
    // 受け皿
    ctx.fillStyle = '#1d3349'; ctx.fillRect(0, wy(sim.h, scale), W, H);
    // 目盛（ノズルからの距離）
    ctx.save();
    ctx.strokeStyle = 'rgba(180,205,230,0.35)'; ctx.fillStyle = 'rgba(190,210,235,0.7)';
    ctx.lineWidth = Math.max(1, scale * 0.0004);
    ctx.font = `${Math.round(scale * 0.0042)}px "IBM Plex Mono", ui-monospace, monospace`;
    const x0 = scale * 0.012;
    for (let cm = 0; cm <= 24; cm++) {
      const y = wy(sim.yN - cm / 100, scale);
      if (y > H) break;
      const len = cm % 5 === 0 ? scale * 0.006 : scale * 0.003;
      ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x0 + len, y); ctx.stroke();
      if (cm % 5 === 0) ctx.fillText(`${cm}`, x0 + len + 3, y + 4);
    }
    ctx.fillText('ノズルから [cm]', x0, wy(sim.yN, scale) - scale * 0.006);
    ctx.restore();

    // サーバー本体とノズル
    const nx = sim.xN * scale, ny = wy(sim.yN, scale);
    ctx.fillStyle = '#d8e0e8';
    ctx.fillRect(nx - scale * 0.03, 0, scale * 0.06, wy(sim.yN + 0.004, scale));
    ctx.fillStyle = '#9aa8b6';
    ctx.fillRect(nx - scale * sim.rNoz, wy(sim.yN + 0.004, scale), 2 * scale * sim.rNoz, scale * 0.004);
    // レバー
    ctx.save();
    ctx.translate(nx + scale * 0.03, wy(sim.yN + 0.012, scale));
    ctx.rotate(sim.lever ? 0.5 : 0);
    ctx.fillStyle = sim.lever ? '#4dabf7' : '#7b8a99';
    ctx.fillRect(0, -scale * 0.003, scale * 0.028, scale * 0.006);
    ctx.restore();
    // ノズル汚染のフラッシュ
    if (now < flashUntil) {
      ctx.strokeStyle = 'rgba(255,90,90,0.95)'; ctx.lineWidth = scale * 0.0015;
      ctx.beginPath(); ctx.arc(nx, ny, scale * 0.012, 0, Math.PI * 2); ctx.stroke();
    }

    // コップ（後ろ側のガラス）
    const c = sim.cup;
    const corners = [[-c.rb, 0], [c.rb, 0], [c.rt, c.hc], [-c.rt, c.hc]].map(([lx, ly]) => c.toWorld(lx, ly));
    ctx.fillStyle = 'rgba(200,225,245,0.07)';
    ctx.beginPath();
    corners.forEach(([x, y], k) => k ? ctx.lineTo(x * scale, wy(y, scale)) : ctx.moveTo(x * scale, wy(y, scale)));
    ctx.closePath(); ctx.fill();

    renderWater(scale);

    // 気泡（見やすさのため実寸の 3 倍で描画）
    ctx.lineWidth = Math.max(0.8, scale * 0.00025);
    for (let k = 0; k < sim.nb; k++) {
      const w = sim.bw[k];
      const a = Math.min(0.9, 0.15 + 0.25 * Math.log10(1 + w));
      if (a < 0.16 && Math.random() > w * 5) continue;
      ctx.strokeStyle = `rgba(235,248,255,${a})`;
      const r = Math.max(1.3, sim.bR[k] * scale * 3);
      ctx.beginPath(); ctx.arc(sim.bx[k] * scale, wy(sim.by[k], scale), r, 0, Math.PI * 2); ctx.stroke();
    }
    // 飛沫（実寸の 6 倍）。代表する実際の数が少ないものは薄く
    for (let k = 0; k < sim.nd; k++) {
      const w = sim.dw[k];
      const a = Math.min(1, 0.12 + 0.3 * Math.log10(1 + w));
      if (w < 0.05 && Math.random() > w * 20) continue;
      ctx.fillStyle = `rgba(255,236,170,${a})`;
      const r = Math.max(1.0, sim.dr[k] * scale * 6);
      const x = sim.dx_[k] * scale, y = wy(sim.dy_[k], scale);
      ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = `rgba(255,236,170,${a * 0.4})`; ctx.lineWidth = r;
      ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x - sim.dvx[k] * scale * 0.004, y + sim.dvy[k] * scale * 0.004); ctx.stroke();
    }

    // コップ（壁）
    ctx.strokeStyle = 'rgba(215,235,250,0.85)'; ctx.lineWidth = c.tw * scale; ctx.lineCap = 'round';
    for (const [ax, ay, bx, by] of c.segs) {
      ctx.beginPath(); ctx.moveTo(ax * scale, wy(ay, scale)); ctx.lineTo(bx * scale, wy(by, scale)); ctx.stroke();
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.35)'; ctx.lineWidth = Math.max(1, c.tw * scale * 0.3);
    const [hx0, hy0] = c.toWorld(c.rb + c.tw * 0.3, 0.01), [hx1, hy1] = c.toWorld(c.rt + c.tw * 0.3, c.hc * 0.9);
    ctx.beginPath(); ctx.moveTo(hx0 * scale, wy(hy0, scale)); ctx.lineTo(hx1 * scale, wy(hy1, scale)); ctx.stroke();

    // 落下距離 H の寸法線
    const st = sim.stats;
    if (sim.lever && st.hasImpact && st.H > 0.005) {
      const yb = sim.yN - st.H;
      const xl = nx - scale * 0.018;
      ctx.strokeStyle = 'rgba(255,214,102,0.8)'; ctx.fillStyle = 'rgba(255,214,102,0.95)';
      ctx.lineWidth = Math.max(1, scale * 0.0004);
      ctx.beginPath(); ctx.moveTo(xl, ny); ctx.lineTo(xl, wy(yb, scale)); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(xl - 4, ny); ctx.lineTo(xl + 4, ny); ctx.moveTo(xl - 4, wy(yb, scale)); ctx.lineTo(xl + 4, wy(yb, scale)); ctx.stroke();
      ctx.font = `${Math.round(scale * 0.0048)}px "IBM Plex Mono", ui-monospace, monospace`;
      ctx.textAlign = 'right';
      ctx.fillText(`H ${(st.H * 100).toFixed(1)} cm`, xl - 6, (ny + wy(yb, scale)) / 2);
      ctx.textAlign = 'left';
    }
    // 再生速度
    ctx.fillStyle = 'rgba(190,210,235,0.75)';
    ctx.font = `${Math.round(scale * 0.0045)}px "IBM Plex Mono", ui-monospace, monospace`;
    ctx.textAlign = 'right';
    ctx.fillText(`${ui.speed === 1 ? '実時間' : 'スロー ' + (1 / ui.speed).toFixed(0) + '倍'}  t = ${sim.time.toFixed(2)} s`, W - scale * 0.006, H - scale * 0.006);
    ctx.textAlign = 'left';
  }

  // ------------------------------------------------------------------ 計測値
  const fmtRate = (x) => x < 1 ? x.toFixed(1) : x < 1e4 ? Math.round(x).toLocaleString('ja-JP') : (x / 1e4).toFixed(x < 1e5 ? 1 : 0) + '万';
  function updateReadouts() {
    const s = sim.stats;
    const active = sim.lever && s.hasImpact;
    $('r-H').textContent = active ? (s.H * 100).toFixed(1) : '—';
    $('r-U').textContent = active ? s.U.toFixed(2) : '—';
    $('r-Ueff').textContent = active ? s.Ueff.toFixed(2) : '—';
    $('r-Ve').textContent = active ? s.Ve.toFixed(2) : '—';
    $('r-mode').textContent = !active ? '注いでいません' : s.wallFrac > 0.5 ? '壁を伝って入る' : '水面に直撃';
    const ratio = active ? s.Ueff / Math.max(s.Ve, 1e-3) : 0;
    const bar = $('gauge-fill');
    bar.style.width = Math.min(100, ratio * 66.7) + '%';
    bar.dataset.state = ratio > 1 ? 'bad' : ratio > 0.85 ? 'warn' : 'ok';
    $('gauge-label').textContent = !active ? '' : ratio > 1 ? '空気を巻き込んでいる' : ratio > 0.85 ? 'ぎりぎり' : '巻き込みなし';
    $('r-bub').textContent = fmtRate(s.bubbleRate);
    $('r-drop').textContent = fmtRate(s.dropRate);
    $('r-noz').textContent = fmtRate(s.nozzleHits + s.bigNozzleHits);
    $('r-rim').textContent = fmtRate(s.rimEscape);
    if (s.nozzleHits + s.bigNozzleHits > lastHits + 0.5) { flashUntil = performance.now() + 400; lastHits = s.nozzleHits + s.bigNozzleHits; }
  }

  // ------------------------------------------------------------------ ループ
  let last = performance.now(), acc = 0, frames = 0, tFps = last;
  function loop(now) {
    const realDt = Math.min((now - last) / 1000, 1 / 20);
    last = now;
    if (ui.playing) {
      const dtSim = realDt * ui.speed;
      moveCupToward(dtSim);
      if (pendingPour > 0) { pendingPour -= dtSim; if (pendingPour <= 0) { sim.resetStats(); lastHits = 0; sim.lever = true; updatePourBtn(); } }
      sim.step(dtSim);
    }
    draw(now);
    updateReadouts();
    frames++;
    if (now - tFps > 1000) { $('fps').textContent = `${frames} fps · 粒子 ${sim.np.toLocaleString('ja-JP')}`; frames = 0; tFps = now; }
    requestAnimationFrame(loop);
  }

  // 起動
  newSim();
  resize();
  syncLabels();
  updatePourBtn();
  window.addEventListener('resize', resize);
  if (window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
    ui.playing = false; $('play').textContent = '再生';
  }
  requestAnimationFrame((t) => { last = t; loop(t); });
  window.__pourApp = { get sim() { return sim; }, ui, applyPreset };
})();
