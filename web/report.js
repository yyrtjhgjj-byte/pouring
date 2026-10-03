/*
 * report.js — レポートのグラフ（SVG・テーマ追従・ホバーで数値表示・表で見る）
 * データは data.js（window.POUR_DATA）から読む。
 */
(function () {
  'use strict';
  const D = window.POUR_DATA;
  if (!D) return;
  const NS = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs = {}, parent) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (parent) parent.appendChild(e);
    return e;
  };
  const fmt = (v, d = 2) => (Math.abs(v) >= 1e4 ? v.toExponential(1).replace('e+', 'e') : v.toFixed(d));

  function niceTicks(lo, hi, n = 5) {
    const span = hi - lo, step0 = span / n;
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || mag * 10;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }

  /**
   * opts: { series: [{name, color, points:[[x,y]], label?, dash?, marker?}], xLabel, yLabel,
   *         xDomain, yDomain, logY, logX, refX:[{x,label}], refY:[{y,label}], band:{x0,x1,label},
   *         diag:{label}, height, xFmt, yFmt, unitX, unitY }
   */
  function lineChart(host, opts) {
    const wrap = document.createElement('div');
    wrap.className = 'chart';
    host.appendChild(wrap);
    if (opts.series.length >= 2) {
      const lg = document.createElement('div');
      lg.className = 'legend';
      opts.series.forEach(s => {
        const it = document.createElement('span');
        it.innerHTML = `<i style="background:var(${s.color})" class="${s.marker ? 'dot' : ''}"></i>${s.name}`;
        lg.appendChild(it);
      });
      wrap.appendChild(lg);
    }
    const box = document.createElement('div');
    box.className = 'plot';
    wrap.appendChild(box);
    const tip = document.createElement('div');
    tip.className = 'tip'; tip.hidden = true;
    box.appendChild(tip);

    const draw = () => {
      box.querySelectorAll('svg').forEach(s => s.remove());
      const W = Math.max(box.clientWidth, 280), H = opts.height || 260;
      const m = { l: 52, r: opts.rightPad ?? 96, t: 12, b: 42 };
      const svg = el('svg', { width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': opts.aria || opts.yLabel });
      box.insertBefore(svg, tip);
      const [x0, x1] = opts.xDomain, [y0, y1] = opts.yDomain;
      const lx = (v) => opts.logX ? Math.log10(v) : v, ly = (v) => opts.logY ? Math.log10(Math.max(v, 1e-12)) : v;
      const X = (v) => m.l + (lx(v) - lx(x0)) / (lx(x1) - lx(x0)) * (W - m.l - m.r);
      const Y = (v) => H - m.b - (ly(v) - ly(y0)) / (ly(y1) - ly(y0)) * (H - m.t - m.b);
      const clipY = (v) => Math.min(Math.max(v, y0), y1);
      // 帯
      if (opts.band) {
        const bx0 = X(Math.max(opts.band.x0, x0)), bx1 = X(Math.min(opts.band.x1, x1));
        el('rect', { x: bx0, y: m.t, width: Math.max(bx1 - bx0, 0), height: H - m.t - m.b, fill: 'var(--band)' }, svg);
        const t = el('text', { x: bx0 + 6, y: m.t + 14, class: 'ann' }, svg); t.textContent = opts.band.label;
      }
      // グリッドと目盛
      const yt = opts.logY ? (() => { const a = []; for (let e = Math.ceil(Math.log10(y0)); e <= Math.log10(y1); e++) a.push(Math.pow(10, e)); return a; })() : (opts.yTicks || niceTicks(y0, y1));
      yt.forEach(v => {
        el('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: 'grid' }, svg);
        const t = el('text', { x: m.l - 6, y: Y(v) + 4, class: 'tick', 'text-anchor': 'end' }, svg);
        t.textContent = opts.yFmt ? opts.yFmt(v) : (opts.logY ? (v >= 1e4 ? `10${sup(Math.log10(v))}` : v.toString()) : v);
      });
      const xt = opts.logX ? (opts.xTicks || [0.02, 0.05, 0.1, 0.2, 0.5, 1]) : (opts.xTicks || niceTicks(x0, x1, 6));
      xt.forEach(v => {
        if (v < x0 || v > x1) return;
        const t = el('text', { x: X(v), y: H - m.b + 16, class: 'tick', 'text-anchor': 'middle' }, svg);
        t.textContent = opts.xFmt ? opts.xFmt(v) : v;
      });
      el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'axis' }, svg);
      const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 6, class: 'axlab', 'text-anchor': 'middle' }, svg); xl.textContent = opts.xLabel;
      const yl = el('text', { x: 12, y: m.t + (H - m.t - m.b) / 2, class: 'axlab', 'text-anchor': 'middle', transform: `rotate(-90 12 ${m.t + (H - m.t - m.b) / 2})` }, svg); yl.textContent = opts.yLabel;
      // 参照線
      (opts.refX || []).forEach(r => {
        el('line', { x1: X(r.x), x2: X(r.x), y1: m.t, y2: H - m.b, class: 'ref' }, svg);
        const t = el('text', { x: X(r.x) + 5, y: m.t + 28, class: 'ann' }, svg); t.textContent = r.label;
      });
      (opts.refY || []).forEach(r => {
        el('line', { x1: m.l, x2: W - m.r, y1: Y(r.y), y2: Y(r.y), class: 'ref' }, svg);
        const t = el('text', { x: W - m.r + 4, y: Y(r.y) + 4, class: 'ann' }, svg); t.textContent = r.label;
      });
      if (opts.diag) {
        const xa = Math.max(x0, y0), xb = Math.min(x1, y1);
        el('line', { x1: X(xa), y1: Y(xa), x2: X(xb), y2: Y(xb), class: 'ref' }, svg);
        const t = el('text', { x: X(xb) + 4, y: Y(xb) + 4, class: 'ann' }, svg); t.textContent = opts.diag.label;
      }
      // 系列
      opts.series.forEach(s => {
        const pts = s.points.filter(p => p[0] >= x0 && p[0] <= x1 && (!opts.logY || p[1] > 0));
        if (!pts.length) return;
        if (!s.marker) {
          const d = pts.map((p, k) => `${k ? 'L' : 'M'}${X(p[0]).toFixed(1)},${Y(clipY(p[1])).toFixed(1)}`).join('');
          el('path', { d, fill: 'none', stroke: `var(${s.color})`, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }, svg);
        } else {
          pts.forEach(p => el('circle', { cx: X(p[0]), cy: Y(clipY(p[1])), r: 5, fill: `var(${s.color})`, stroke: 'var(--panel)', 'stroke-width': 2 }, svg));
        }
        if (s.label !== false) {
          const p = pts[pts.length - 1];
          const t = el('text', { x: X(p[0]) + 8, y: Y(clipY(p[1])) + 4 + (s.labelDy || 0), class: 'dlab' }, svg);
          t.textContent = s.label || s.name;
        }
      });
      // ホバー
      const hit = el('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, svg);
      const cross = el('line', { y1: m.t, y2: H - m.b, class: 'cross', visibility: 'hidden' }, svg);
      const move = (ev) => {
        const r = svg.getBoundingClientRect();
        const px = ev.clientX - r.left;
        const fx = (px - m.l) / (W - m.l - m.r);
        const xv = opts.logX ? Math.pow(10, lx(x0) + fx * (lx(x1) - lx(x0))) : x0 + fx * (x1 - x0);
        let html = `<b>${opts.xFmt ? opts.xFmt(xv) : xv.toFixed(1)}${opts.unitX || ''}</b>`;
        let near = null;
        opts.series.forEach(s => {
          let best = null, bd = Infinity;
          s.points.forEach(p => { const d = Math.abs(lx(p[0]) - lx(xv)); if (d < bd) { bd = d; best = p; } });
          if (best) { html += `<br><i style="background:var(${s.color})"></i>${s.name}: ${opts.yFmt ? opts.yFmt(best[1]) : fmt(best[1], opts.yDigits ?? 2)}${opts.unitY || ''}`; near = best; }
        });
        cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('visibility', 'visible');
        tip.innerHTML = html; tip.hidden = false;
        const tw = tip.offsetWidth;
        tip.style.left = Math.min(Math.max(px + 12, 0), W - tw - 4) + 'px';
        tip.style.top = (m.t + 4) + 'px';
      };
      hit.addEventListener('pointermove', move);
      hit.addEventListener('pointerleave', () => { tip.hidden = true; cross.setAttribute('visibility', 'hidden'); });
    };
    draw();
    let rt;
    new ResizeObserver(() => { clearTimeout(rt); rt = setTimeout(draw, 80); }).observe(box);
    // 表で見る
    if (opts.table) {
      const det = document.createElement('details');
      det.className = 'tableview';
      det.innerHTML = `<summary>表で見る</summary>`;
      const sc = document.createElement('div'); sc.className = 'tscroll';
      const tb = document.createElement('table');
      const hdr = `<tr><th>${opts.xLabel}</th>${opts.series.map(s => `<th>${s.name}</th>`).join('')}</tr>`;
      const xs = opts.table.xs;
      const rows = xs.map(x => `<tr><td>${opts.xFmt ? opts.xFmt(x) : x}</td>${opts.series.map(s => {
        const p = s.points.find(q => Math.abs(q[0] - x) < 1e-9);
        return `<td>${p ? (opts.yFmt ? opts.yFmt(p[1]) : fmt(p[1], opts.yDigits ?? 2)) : '—'}</td>`;
      }).join('')}</tr>`).join('');
      tb.innerHTML = hdr + rows;
      sc.appendChild(tb); det.appendChild(sc); wrap.appendChild(det);
    }
    return wrap;
  }
  function sup(n) { return String(n).split('').map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+c] || c).join(''); }
  const human = (v) => v >= 1e8 ? (v / 1e8).toFixed(1) + '億' : v >= 1e4 ? (v / 1e4).toFixed(v >= 1e5 ? 0 : 1) + '万' : v >= 10 ? Math.round(v).toString() : v >= 1 ? v.toFixed(1) : v.toFixed(2);

  // ------------------------------------------------------------------ 1. 着水速度と巻き込み開始速度
  const host1 = document.getElementById('chart-cross');
  if (host1) {
    const P = D.model.plunge;
    const Hs = D.model.H_star * 100;
    lineChart(host1, {
      series: [
        { name: '水柱が水面に当たる速さ U', color: '--s-orange', points: P.map(r => [r.H * 100, r.U]), label: '当たる速さ' },
        { name: '空気を巻き込み始める速さ Vₑ', color: '--s-aqua', points: P.map(r => [r.H * 100, r.Ve]), label: '巻き込み開始' },
      ],
      xDomain: [2, 30], yDomain: [0, 3], xLabel: 'ノズルから水面までの落下距離 H [cm]', yLabel: '速さ [m/s]',
      band: { x0: Hs, x1: 30, label: '空気を巻き込む' }, refX: [{ x: Hs, label: `H* ≈ ${Hs.toFixed(0)} cm` }],
      unitX: ' cm', unitY: ' m/s', xFmt: v => v.toFixed(0), aria: '落下距離と速さのグラフ',
      table: { xs: [2, 4, 6, 8, 10, 12, 15, 20, 25, 30] },
    });
  }

  // ------------------------------------------------------------------ 2. 気泡の数（直撃 vs 壁沿い）
  const host2 = document.getElementById('chart-bubbles');
  if (host2) {
    const P = D.model.plunge, W30 = D.model.wall['30'], W4 = D.model.wall['4'];
    lineChart(host2, {
      series: [
        { name: 'まっすぐ水面に落とす', color: '--s-orange', points: P.map(r => [r.H * 100, Math.max(r.bubble_rate, 0.1)]), label: '水面に直撃' },
        { name: '直立のまま内壁をかすめる', color: '--s-aqua', points: W4.map(r => [r.H * 100, Math.max(r.bubble_rate, 0.1)]), label: '内壁かすめ', labelDy: -8 },
        { name: '30°傾けて内壁に当てる', color: '--s-blue', points: W30.map(r => [r.H * 100, Math.max(r.bubble_rate, 0.1)]), label: '30°傾け', labelDy: 8 },
      ],
      xDomain: [2, 30], yDomain: [0.1, 1e7], logY: true,
      xLabel: '落下距離 H [cm]（壁沿いは壁に当たるまで）', yLabel: 'できる気泡 [個/秒]',
      yFmt: human, unitX: ' cm', unitY: ' 個/秒', xFmt: v => v.toFixed(0), aria: '落下距離と気泡の数',
      table: { xs: [4, 6, 8, 10, 12, 15, 20, 25, 30] },
    });
  }

  // ------------------------------------------------------------------ 3. 水量ごとの傾き
  const host3 = document.getElementById('chart-schedule');
  if (host3) {
    const S = D.model.schedule.filter(r => r.theta_max !== null);
    lineChart(host3, {
      series: [
        { name: 'こぼさず内壁を狙える上限', color: '--s-blue', points: S.map(r => [r.V, r.theta_max]), label: '上限' },
        { name: 'おすすめ（上限 − 5°）', color: '--s-orange', points: S.map(r => [r.V, r.theta_rec]), label: 'おすすめ' },
      ],
      xDomain: [40, 300], yDomain: [0, 50], xLabel: `コップの中の水 [mL]（${D.model.cup.volume_mL.toFixed(0)} mL のタンブラー）`, yLabel: '傾き [°]',
      band: { x0: S[S.length - 1].V + 5, x1: 300, label: '直立でよい（落下距離が短い）' },
      unitX: ' mL', unitY: '°', xFmt: v => v.toFixed(0), yDigits: 0, aria: '水量と傾き',
      table: { xs: S.map(r => r.V).filter((v, i) => i % 2 === 0) },
    });
  }

  // ------------------------------------------------------------------ 4. 気泡の大きさ → 飛沫の到達高さ
  const host4 = document.getElementById('chart-jetdrops');
  if (host4) {
    const J = D.model.jetdrops.filter(r => r.n > 0);
    lineChart(host4, {
      series: [
        { name: '空気抵抗あり（実際）', color: '--s-blue', points: J.map(r => [r.R_um / 1000, r.h_cm]), label: '空気抵抗あり' },
      ],
      xDomain: [0.02, 1.6], yDomain: [0, 8], logX: true, xTicks: [0.02, 0.05, 0.1, 0.2, 0.5, 1],
      xLabel: '弾ける気泡の半径 [mm]', yLabel: '一番上の飛沫が届く高さ [cm]',
      xFmt: v => v < 0.1 ? v.toFixed(2) : v.toFixed(1), unitX: ' mm', unitY: ' cm', yDigits: 1, aria: '気泡の大きさと飛沫の高さ',
      table: { xs: J.filter((r, i) => i % 4 === 0).map(r => r.R_um / 1000) },
    });
  }

  // ------------------------------------------------------------------ 5. DNS: 落下距離と水滴の最高到達点
  const host5 = document.getElementById('chart-dns');
  if (host5 && D.dns && D.dns.startup && D.dns.startup.length) {
    const R = D.dns.startup.slice().sort((a, b) => a.H - b.H);
    lineChart(host5, {
      series: [
        { name: '注ぎ始めに打ち上がる水滴の最高到達点（DNS）', color: '--s-orange', points: R.map(r => [r.H * 100, r.apex * 100]), marker: true, label: false },
      ],
      xDomain: [0, 22], yDomain: [0, Math.max(22, ...R.map(r => r.apex * 100 + 2))],
      diag: { label: 'ノズルの高さ' }, rightPad: 90,
      xLabel: 'ノズルから水面までの落下距離 H [cm]', yLabel: '水面からの高さ [cm]',
      unitX: ' cm', unitY: ' cm', xFmt: v => v.toFixed(0), yDigits: 1, aria: 'DNS の水滴到達高さ',
      table: { xs: R.map(r => r.H * 100) },
    });
  }
  const host6 = document.getElementById('chart-dns-air');
  if (host6 && D.dns && D.dns.startup && D.dns.startup.length) {
    const R = D.dns.startup.slice().sort((a, b) => a.H - b.H);
    lineChart(host6, {
      series: [
        { name: '注ぎ始めに閉じ込められる空気（DNS）', color: '--s-blue', points: R.map(r => [r.H * 100, r.air_mL]), marker: true, label: false },
      ],
      xDomain: [0, 22], yDomain: [0, Math.max(1, ...R.map(r => r.air_mL)) * 1.15],
      xLabel: 'ノズルから水面までの落下距離 H [cm]', yLabel: '閉じ込められた空気 [mL]',
      unitX: ' cm', unitY: ' mL', xFmt: v => v.toFixed(0), yDigits: 2, aria: 'DNS の巻き込み空気量', rightPad: 24,
      table: { xs: R.map(r => r.H * 100) },
    });
  }
})();
