// Node でプリセットを走らせて統計を比較する
const { PourSim } = require('../sim.js');
const presets = {
  straight_far: { gap: 0.12, th: 0, aim: 'center' },
  straight_mid: { gap: 0.04, th: 0, aim: 'center' },
  straight_near: { gap: 0.008, th: 0, aim: 'center' },
  tilt_wall: { gap: 0.025, th: 30, aim: 'wall' },
};
const which = process.argv.slice(2);
for (const [name, ps] of Object.entries(presets)) {
  if (which.length && !which.includes(name)) continue;
  const sim = new PourSim({ h: 0.00125, fillFrac: 0.3 });
  const pl = sim.targetPose(ps.gap, ps.th, ps.aim, 0.012);
  sim.teleportCup(pl.x, pl.y, pl.th);
  for (let f = 0; f < 150; f++) sim.step(1 / 240);
  sim.resetStats(); sim.lever = true;
  const t0 = Date.now(); let steps = 0;
  const T = 2.0, fdt = 1 / 240;
  for (let f = 0; f * fdt < T; f++) steps += sim.step(fdt);
  const s = sim.stats;
  console.log(`${name.padEnd(14)} tilt=${pl.deg}° H=${(s.H*100).toFixed(1)}cm U=${s.U.toFixed(2)} Ueff=${s.Ueff.toFixed(2)} Ve=${s.Ve.toFixed(2)} air=${s.air.toExponential(1)} wall=${s.wallFrac.toFixed(2)} bubbles=${s.bubbles.toExponential(1)} drops=${s.drops.toExponential(1)} nozzle=${s.nozzleHits.toFixed(2)} rim=${s.rimEscape.toExponential(1)} bigNoz=${s.bigNozzleHits} spilled=${(s.spilled*1e6).toFixed(1)}mL  (${((Date.now()-t0)/steps).toFixed(1)} ms/substep)`);
}
