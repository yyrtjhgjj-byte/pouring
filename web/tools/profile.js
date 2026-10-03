// 1 フレームあたりの各処理の時間を測る（開発用）
const { PourSim } = require('../sim.js');
let t0 = Date.now();
const sim = new PourSim({ h: +(process.argv[2] || 0.00125) });
console.log('construct', Date.now() - t0, 'ms  np', sim.np, 'rest', sim.restDensity.toFixed(2), 'grid', sim.nx, 'x', sim.ny);
const pl = sim.targetPose(0.12, 0, 'center', 0.02);
sim.teleportCup(pl.x, pl.y, pl.th);
const fns = ['emit', 'integrate', 'pushApart', 'collide', 'transfer', 'updateDensity', 'solvePressure', 'detectImpacts', 'stepBubbles', 'stepDrops'];
const T = {};
for (const fn of fns) {
  const orig = sim[fn].bind(sim); T[fn] = 0;
  sim[fn] = (...a) => { const t = process.hrtime.bigint(); const r = orig(...a); T[fn] += Number(process.hrtime.bigint() - t) / 1e6; return r; };
}
sim.lever = true;
t0 = Date.now();
let n = 0;
for (let f = 0; f < 20; f++) n += sim.step(1 / 240);
console.log('20 frames', n, 'substeps', Date.now() - t0, 'ms  np', sim.np);
for (const fn of fns) console.log(fn.padEnd(14), T[fn].toFixed(1), 'ms');
