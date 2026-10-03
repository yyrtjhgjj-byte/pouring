// プリセットを走らせて粒子位置を JSON に書き出す（matplotlib で確認用）
const { PourSim } = require('../sim.js');
const fs = require('fs');
const [gap, th, aim, T, out, fill] = [+process.argv[2], +process.argv[3], process.argv[4], +process.argv[5], process.argv[6], +(process.argv[7] || 0.3)];
const sim = new PourSim({ h: 0.00125, fillFrac: fill });
const pl = sim.targetPose(gap, th, aim, 0.012);
console.log('pose', pl);
sim.teleportCup(pl.x, pl.y, pl.th);
for (let f = 0; f < 150; f++) sim.step(1 / 240);
sim.resetStats(); sim.lever = true;
const frames = [];
const t0 = Date.now(); let ns = 0;
for (let f = 0; f * (1 / 240) < T; f++) {
  ns += sim.step(1 / 240);
  if (f % 36 === 0) {
    const P = [];
    for (let i = 0; i < sim.np; i++) P.push([+sim.pos[2*i].toFixed(5), +sim.pos[2*i+1].toFixed(5), sim.state[i]]);
    const B = []; for (let k = 0; k < sim.nb; k++) B.push([sim.bx[k], sim.by[k]]);
    const D = []; for (let k = 0; k < sim.nd; k++) D.push([sim.dx_[k], sim.dy_[k]]);
    frames.push({ t: sim.time, P, B, D, segs: sim.cup.segs, stats: { ...sim.stats, events: undefined } });
  }
}
console.log('ms/substep', ((Date.now()-t0)/ns).toFixed(2), 'final', JSON.stringify({ ...sim.stats, events: undefined }));
fs.writeFileSync(out, JSON.stringify({ xN: sim.xN, yN: sim.yN, W: sim.W, H: sim.H, frames }));
