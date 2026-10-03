// ヘッドレス Chromium でページを開き、スクリーンショットとコンソールエラーを確認する
const { chromium } = require(process.env.PW || 'playwright');
(async () => {
  const url = process.argv[2], out = process.argv[3], wait = +(process.argv[4] || 4000), preset = process.argv[5];
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: +(process.env.VW || 1280), height: +(process.env.VH || 900) }, deviceScaleFactor: 1 });
  const errs = [];
  page.on('console', m => { if (m.type() === 'error' || m.type() === 'warning') errs.push(m.type() + ': ' + m.text()); });
  page.on('pageerror', e => errs.push('pageerror: ' + e.message));
  await page.goto(url);
  if (process.env.SCHEME) await page.emulateMedia({ colorScheme: process.env.SCHEME });
  await page.waitForTimeout(500);
  if (preset) await page.click(`.preset[data-preset="${preset}"]`);
  await page.waitForTimeout(wait);
  const info = await page.evaluate(() => { const a = window.__pourApp; if (!a) return null; const s = a.sim.stats; return { t: a.sim.time, np: a.sim.np, H: s.H, Ueff: s.Ueff, Ve: s.Ve, bub: s.bubbleRate, drop: s.dropRate, fps: document.getElementById('fps').textContent }; });
  console.log(JSON.stringify(info));
  await page.screenshot({ path: out, fullPage: true });
  console.log(errs.join('\n') || 'no console errors');
  await browser.close();
})();
