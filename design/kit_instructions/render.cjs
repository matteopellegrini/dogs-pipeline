const puppeteer = require('/Users/matteopellegrini/Downloads/dogs/dogs-app/node_modules/puppeteer-core');
(async () => {
  const b = await puppeteer.launch({ executablePath: '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', headless: 'new' });
  const p = await b.newPage();
  await p.goto('file://' + process.env.S + '/sheet.html', { waitUntil: 'networkidle0' });
  await p.pdf({ path: process.env.S + '/ProsperK9_Kit_Instructions_v4.pdf', format: 'Letter', printBackground: true, preferCSSPageSize: true });
  await p.setViewport({ width: 816, height: 1056, deviceScaleFactor: 2 });
  await p.screenshot({ path: process.env.S + '/sheet_preview.png', fullPage: false });
  await b.close();
})().catch(e => { console.error('ERR', e.message); process.exit(1); });
