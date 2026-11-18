// Numerical Techniques Labs showcase: every lab in sripts/ is compiled unchanged to WebAssembly
// (web/build.sh); worker.js runs its main() and this page parses what it printed into tables
// and charts. Curves badged "JS" are reference values computed here, not by the labs.

const REPO = 'https://github.com/Andrew1407/NumericalTechniques';
const SRC = `${REPO}/blob/main/sripts`;
const $ = (s, root = document) => root.querySelector(s);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (v, d = 4) => {
  if (typeof v !== 'number') return String(v);
  if (Number.isNaN(v)) return 'nan';
  if (!Number.isFinite(v)) return v > 0 ? 'inf' : '-inf';
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(2);
  return String(+v.toFixed(d));
};
const range = (n, f = (i) => i) => Array.from({ length: n }, (_, i) => f(i));
const linspace = (a, b, n) => range(n, (i) => a + ((b - a) * i) / (n - 1));
const stat = (value, label, cls = '') => `<div class="stat ${cls}"><b>${value}</b><span>${label}</span></div>`;
const JS = '<span class="lang ref" title="Reference value computed by this page in JavaScript, not by the lab">JS</span>';
// C printf spellings of special values: inf, -inf, nan, -nan
const num = (s) => {
  const t = String(s).trim().replace(/[;:,]$/, '');
  if (/^-?nan$/i.test(t)) return NaN;
  if (/^inf$/i.test(t)) return Infinity;
  if (/^-inf$/i.test(t)) return -Infinity;
  return Number(t);
};

// ---------------------------------------------------------------- runtime --
// One fresh worker (and wasm instance) per run, so a lab that never stops can be cut off.
const TIME_LIMIT = 8000;
function runWasm(lab, globals) {
  return new Promise((resolve, reject) => {
    const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    const timer = setTimeout(() => { w.terminate(); reject(new Error(`The program didn't finish within ${TIME_LIMIT / 1000} s, so it was stopped.`)); }, TIME_LIMIT);
    w.onmessage = ({ data }) => { clearTimeout(timer); w.terminate(); data.ok ? resolve(data) : reject(new Error(data.error)); };
    w.onerror = (e) => { clearTimeout(timer); w.terminate(); reject(new Error(e.message || 'The WebAssembly worker failed to start.')); };
    w.postMessage({ lab, globals });
  });
}

// ------------------------------------------------------------------ charts --
// spec: { title, x, y, series: [{ name, xs, ys, color, dash, dots, line, width, hidden }],
//         height, log, logToggle, view: { y: [lo, hi] }, areas: [{ xs, ys, color }],
//         segments: [{ pts: [[x, y], …], color, dash }], marks: [{ x, y, label, color }] }
const charts = new Map();
let chartSeq = 0;

function niceStep(span, count) {
  const raw = span / Math.max(1, count);
  const p = 10 ** Math.floor(Math.log10(raw || 1));
  const m = raw / p;
  return (m >= 5 ? 10 : m >= 2 ? 5 : m >= 1 ? 2 : 1) * p;
}
function ticks(min, max, count) {
  if (min === max) { min -= Math.abs(min) * 0.1 || 1; max += Math.abs(max) * 0.1 || 1; }
  const step = niceStep(max - min, count);
  const lo = Math.floor(min / step) * step, hi = Math.ceil(max / step) * step;
  const out = [];
  for (let v = lo; v <= hi + step / 2; v += step) out.push(+v.toPrecision(12));
  return { lo, hi, ticks: out };
}

function chart(spec) {
  const id = `c${++chartSeq}`;
  charts.set(id, { spec, hidden: new Set(spec.series.flatMap((s, i) => (s.hidden ? [i] : []))), log: !!spec.log, zoom: null, full: false });
  const keys = spec.series.length > 1
    ? `<span class="keys">${spec.series.map((s, i) => `<button type="button" data-key="${i}" aria-pressed="${!s.hidden}" title="Show or hide this series"><i style="background:var(${s.color})${s.dash ? ';opacity:.7' : ''}"></i>${esc(s.name)}${s.ref ? '&nbsp;' + JS : ''}</button>`).join('')}</span>` : '';
  const tools = (spec.logToggle ? `<button type="button" class="mini" data-log aria-pressed="${!!spec.log}" title="Logarithmic y axis">Log y</button>` : '')
    + (spec.view ? `<button type="button" class="mini" data-full aria-pressed="false" title="Show the full y range, including the values the default view clips">Full range</button>` : '')
    + '<button type="button" class="mini" data-reset hidden title="Back to the default view (or double-click the chart)">Reset zoom</button>';
  return `<figure class="chart" data-chart="${id}" data-copyable data-title="${esc(spec.title)}"><figcaption><span>${esc(spec.title)}</span>${keys}<span class="copy-actions">${tools}${actionButtons('png', 'csv')}</span></figcaption><div class="plot"></div>${spec.note ? `<div class="chart-hint">${spec.note}</div>` : ''}</figure>`;
}

function drawCharts(root = document) {
  root.querySelectorAll('figure.chart').forEach((fig) => { if ($('.plot', fig).offsetParent) drawChart(fig); });
}
let resizeTimer;
window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(() => drawCharts(), 120); });

function drawChart(fig) {
  const id = fig.dataset.chart, st = charts.get(id), spec = st.spec, plot = $('.plot', fig);
  const W = Math.max(260, plot.clientWidth), H = spec.height || 250;
  const M = { l: 60, r: 14, t: 10, b: 38 };
  const log = st.log;
  const T = log ? Math.log10 : (v) => v;
  const okY = (y) => Number.isFinite(y) && (!log || y > 0);
  const vis = spec.series.map((s, i) => ({ ...s, i })).filter((s) => !st.hidden.has(s.i));
  const pts = vis.flatMap((s) => s.xs.map((x, k) => [x, s.ys[k]])).concat((spec.areas || []).flatMap((a) => a.xs.map((x, k) => [x, a.ys[k]])))
    .filter(([x, y]) => Number.isFinite(x) && okY(y));
  $('[data-reset]', fig).hidden = !st.zoom;
  if (!pts.length) { plot.innerHTML = '<p class="hint">No finite data to plot.</p>'; return; }
  let x0, x1, y0, y1, yt;
  const cnt = Math.max(3, Math.floor((W - M.l) / 90));
  if (st.zoom) {
    ({ x0, x1, y0, y1 } = st.zoom);
  } else {
    x0 = Math.min(...pts.map((p) => p[0])); x1 = Math.max(...pts.map((p) => p[0]));
    if (spec.xPad) { const d = (x1 - x0) * spec.xPad; x0 -= d; x1 += d; }
    let ys = pts.map((p) => T(p[1]));
    if (spec.areas && !log) ys = ys.concat(0);
    let lo = Math.min(...ys), hi = Math.max(...ys);
    if (spec.view && !st.full && !log) { lo = Math.max(lo, spec.view.y[0]); hi = Math.min(hi, spec.view.y[1]); }
    if (log) { y0 = Math.floor(lo); y1 = Math.ceil(hi); if (y0 === y1) { y0 -= 1; y1 += 1; } }
    else { const r = ticks(lo, hi, 5); y0 = r.lo; y1 = r.hi; }
  }
  if (log) {
    const span = y1 - y0, every = Math.max(1, Math.ceil(span / 6));
    yt = range(Math.floor(span) + 2, (k) => Math.ceil(y0) + k).filter((e) => e <= y1 + 1e-9 && e % every === 0);
  } else yt = ticks(y0, y1, 5).ticks.filter((v) => v >= y0 - 1e-9 && v <= y1 + 1e-9);
  const X = (v) => M.l + ((v - x0) / (x1 - x0 || 1)) * (W - M.l - M.r);
  const Yt = (t) => H - M.b - ((t - y0) / (y1 - y0 || 1)) * (H - M.t - M.b);
  const Y = (v) => Yt(T(v));
  const xt = ticks(x0, x1, cnt).ticks.filter((v) => v >= x0 - 1e-9 && v <= x1 + 1e-9 && (!spec.intX || Number.isInteger(v)));
  const ylab = (t) => (log ? (Math.abs(t) <= 3 ? fmt(10 ** t, 6) : `1e${t}`) : fmt(t, 6));
  const grid = yt.map((t) => `<line x1="${M.l}" x2="${W - M.r}" y1="${Yt(t)}" y2="${Yt(t)}" stroke="var(--grid)"/><text x="${M.l - 8}" y="${Yt(t)}" text-anchor="end" dominant-baseline="middle">${ylab(t)}</text>`).join('')
    + xt.map((v) => `<line x1="${X(v)}" x2="${X(v)}" y1="${H - M.b}" y2="${H - M.b + 4}" stroke="var(--muted)"/><text x="${X(v)}" y="${H - M.b + 16}" text-anchor="middle">${fmt(v, 6)}</text>`).join('');
  const axes = `<line x1="${M.l}" x2="${W - M.r}" y1="${H - M.b}" y2="${H - M.b}" stroke="var(--line)"/>`
    + (!log && y0 < 0 && y1 > 0 ? `<line x1="${M.l}" x2="${W - M.r}" y1="${Y(0)}" y2="${Y(0)}" stroke="var(--muted)" stroke-opacity=".45"/>` : '')
    + `<text x="${(M.l + W - M.r) / 2}" y="${H - 4}" text-anchor="middle">${esc(spec.x)}</text>`
    + `<text transform="translate(12 ${(M.t + H - M.b) / 2}) rotate(-90)" text-anchor="middle">${esc(spec.y)}${log ? ' (log)' : ''}</text>`;
  const P = (x, y) => `${X(x).toFixed(1)},${Y(y).toFixed(1)}`;
  const areas = log ? '' : (spec.areas || []).map((a) => range(a.xs.length - 1, (k) => `<polygon points="${P(a.xs[k], 0)} ${P(a.xs[k], a.ys[k])} ${P(a.xs[k + 1], a.ys[k + 1])} ${P(a.xs[k + 1], 0)}" fill="var(${a.color})" fill-opacity=".16" stroke="var(${a.color})" stroke-width="1.2"/>`).join('')).join('');
  const segs = (spec.segments || []).map((g) => `<polyline points="${g.pts.filter(([x, y]) => Number.isFinite(x) && okY(y)).map(([x, y]) => P(x, y)).join(' ')}" fill="none" stroke="var(${g.color})" stroke-width="1.2" stroke-opacity=".8" ${g.dash ? 'stroke-dasharray="4 3"' : ''}/>`).join('');
  const lines = vis.map((s) => {
    // break the line at inf / nan (and at y <= 0 on a log axis)
    const runs = [[]];
    s.xs.forEach((x, k) => { if (Number.isFinite(x) && okY(s.ys[k])) runs[runs.length - 1].push(P(x, s.ys[k])); else if (runs[runs.length - 1].length) runs.push([]); });
    const line = s.line === false ? '' : runs.filter((r) => r.length > 1).map((r) => `<polyline points="${r.join(' ')}" fill="none" stroke="var(${s.color})" stroke-width="${s.width || 1.8}" stroke-linejoin="round" ${s.dash ? 'stroke-dasharray="6 4"' : ''}/>`).join('');
    const dots = s.dots || spec.dots || s.line === false || s.xs.length <= 30
      ? s.xs.map((x, k) => (Number.isFinite(x) && okY(s.ys[k]) ? `<circle cx="${X(x)}" cy="${Y(s.ys[k])}" r="${s.r || 2.6}" fill="var(${s.color})"/>` : '')).join('') : '';
    return line + dots;
  }).join('');
  const marks = (spec.marks || []).filter((m) => okY(m.y)).map((m) => `<circle cx="${X(m.x)}" cy="${Y(m.y)}" r="5" fill="var(--panel)" stroke="var(${m.color})" stroke-width="2.4"/><text x="${X(m.x) > W * 0.6 ? X(m.x) - 9 : X(m.x) + 9}" y="${Y(m.y) - 9}" text-anchor="${X(m.x) > W * 0.6 ? 'end' : 'start'}" style="fill:var(--ink)">${esc(m.label)}</text>`).join('');
  const clip = `clip-${id}`;
  plot.innerHTML = `<svg viewBox="0 0 ${W} ${H}" height="${H}" role="img" aria-label="${esc(spec.title)}"><defs><clipPath id="${clip}"><rect x="${M.l - 4}" y="${M.t - 4}" width="${W - M.l - M.r + 8}" height="${H - M.t - M.b + 8}"/></clipPath></defs>${grid}${axes}<g clip-path="url(#${clip})">${areas}${segs}${lines}${marks}</g><line class="cross" y1="${M.t}" y2="${H - M.b}" stroke="var(--muted)" stroke-dasharray="3 3" visibility="hidden"/><g class="marks"></g><rect class="zoom" visibility="hidden" x="0" y="0" width="0" height="0"/></svg><div class="tip" hidden></div>`;

  const svg = $('svg', plot), tip = $('.tip', plot), cross = $('.cross', svg), hl = $('.marks', svg), zr = $('.zoom', svg);
  const toSvg = (e) => { const r = svg.getBoundingClientRect(); return [((e.clientX - r.left) / r.width) * W, ((e.clientY - r.top) / r.height) * H]; };
  const invX = (px) => x0 + ((px - M.l) / (W - M.l - M.r)) * (x1 - x0);
  const invT = (py) => y0 + ((H - M.b - py) / (H - M.t - M.b)) * (y1 - y0);
  const nearest = (xs, v) => {
    let best = -1, bd = Infinity;
    xs.forEach((x, k) => { const d = Math.abs(x - v); if (Number.isFinite(x) && d < bd) { bd = d; best = k; } });
    return best;
  };
  const hide = () => { tip.hidden = true; cross.setAttribute('visibility', 'hidden'); hl.innerHTML = ''; };
  const hover = (e) => {
    const [px] = toSvg(e);
    if (px < M.l - 6 || px > W - M.r + 6 || !vis.length) return hide();
    const v = invX(px);
    const hits = vis.map((s) => { const k = nearest(s.xs, v); return { s, x: s.xs[k], y: s.ys[k] }; }).filter((h) => h.x !== undefined);
    if (!hits.length) return hide();
    hits.sort((a, b) => Math.abs(a.x - v) - Math.abs(b.x - v));
    const hx = X(hits[0].x);
    const shown = hits.filter((h) => Math.abs(X(h.x) - hx) < 1.5);
    cross.setAttribute('x1', hx); cross.setAttribute('x2', hx); cross.setAttribute('visibility', 'visible');
    hl.innerHTML = shown.filter((h) => okY(h.y) && T(h.y) >= y0 && T(h.y) <= y1).map((h) => `<circle cx="${X(h.x)}" cy="${Y(h.y)}" r="4" fill="var(--panel)" stroke="var(${h.s.color})" stroke-width="2"/>`).join('');
    tip.innerHTML = `<div class="tx">${esc(spec.x)} = ${fmt(hits[0].x, 6)}</div>` + shown.map((h) => `<div><i style="background:var(${h.s.color})"></i>${esc(h.s.name || spec.y)}: ${fmt(h.y, 7)}</div>`).join('');
    tip.hidden = false;
    const rect = svg.getBoundingClientRect(), left = (hx / W) * rect.width, tw = tip.offsetWidth;
    tip.style.left = `${left + 12 + tw > rect.width ? Math.max(0, left - tw - 12) : left + 12}px`;
    tip.style.top = `${M.t}px`;
  };
  // drag a box to zoom (mouse and pen); touch keeps scrolling the page and shows the tooltip
  let drag = null;
  svg.addEventListener('pointerdown', (e) => {
    hover(e);
    if (e.pointerType === 'touch' || e.button !== 0) return;
    drag = toSvg(e); svg.setPointerCapture(e.pointerId);
  });
  svg.addEventListener('pointermove', (e) => {
    if (!drag) return hover(e);
    const [px, py] = toSvg(e);
    const clampX = (v) => Math.min(W - M.r, Math.max(M.l, v)), clampY = (v) => Math.min(H - M.b, Math.max(M.t, v));
    const ax = clampX(drag[0]), bx = clampX(px), ay = clampY(drag[1]), by = clampY(py);
    zr.setAttribute('x', Math.min(ax, bx)); zr.setAttribute('width', Math.abs(bx - ax));
    zr.setAttribute('y', Math.min(ay, by)); zr.setAttribute('height', Math.abs(by - ay));
    zr.setAttribute('visibility', Math.abs(bx - ax) > 6 && Math.abs(by - ay) > 6 ? 'visible' : 'hidden');
    hover(e);
  });
  svg.addEventListener('pointerup', (e) => {
    if (!drag) return;
    const vis = zr.getAttribute('visibility') === 'visible';
    drag = null;
    if (!vis) return;
    const zx = +zr.getAttribute('x'), zw = +zr.getAttribute('width'), zy = +zr.getAttribute('y'), zh = +zr.getAttribute('height');
    st.zoom = { x0: invX(zx), x1: invX(zx + zw), y0: invT(zy + zh), y1: invT(zy) };
    drawChart(fig);
  });
  svg.addEventListener('pointerleave', () => { if (!drag) hide(); });
  svg.addEventListener('dblclick', () => { if (st.zoom) { st.zoom = null; drawChart(fig); } });
}

document.addEventListener('click', (e) => {
  const fig = e.target.closest('figure.chart');
  if (!fig) return;
  const st = charts.get(fig.dataset.chart);
  const key = e.target.closest('.keys button[data-key]');
  if (key) {
    const i = +key.dataset.key;
    st.hidden.has(i) ? st.hidden.delete(i) : st.hidden.add(i);
    key.setAttribute('aria-pressed', !st.hidden.has(i));
    return drawChart(fig);
  }
  const b = e.target.closest('button[data-log], button[data-full], button[data-reset]');
  if (!b) return;
  if (b.hasAttribute('data-log')) { st.log = !st.log; st.zoom = null; b.setAttribute('aria-pressed', st.log); }
  if (b.hasAttribute('data-full')) { st.full = !st.full; st.zoom = null; b.setAttribute('aria-pressed', st.full); }
  if (b.hasAttribute('data-reset')) st.zoom = null;
  drawChart(fig);
});

// ------------------------------------------------------------------ tables --
function table(title, head, rows, { cls = '', note = '' } = {}) {
  return `<div data-copyable data-title="${esc(title)}"><div class="section-h"><h2>${esc(title)}${note ? ` <small>· ${note}</small>` : ''}</h2>${actions('table')}</div>
    <div class="table-wrap"><table class="${cls}"><thead><tr>${head.map((h) => `<th class="num">${h}</th>`).join('')}</tr></thead><tbody>${rows.map((r) => `<tr>${r.map((c) => `<td class="num">${typeof c === 'number' ? fmt(c, 7) : c}</td>`).join('')}</tr>`).join('')}</tbody></table></div></div>`;
}

// ------------------------------------------------------------------ maths --
function gauss(A, b) {
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    if (Math.abs(M[c][c]) < 1e-14) return null;
    for (let r = c + 1; r < n; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
  }
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) { let v = M[i][n]; for (let k = i + 1; k < n; k++) v -= M[i][k] * x[k]; x[i] = v / M[i][i]; }
  return x;
}
// least squares fit y ≈ Σ c_j φ_j(x)
function leastSquares(xs, ys, basis) {
  const A = basis.map((f) => basis.map((g) => xs.reduce((s, x) => s + f(x) * g(x), 0)));
  const b = basis.map((f) => xs.reduce((s, x, i) => s + f(x) * ys[i], 0));
  return gauss(A, b);
}
const simpson = (f, a, b, n = 2000) => { const h = (b - a) / n; let s = f(a) + f(b); for (let i = 1; i < n; i++) s += f(a + i * h) * (i % 2 ? 4 : 2); return (s * h) / 3; };

// ---------------------------------------------------------------- parsing --
const kv = (line) => { const m = /^\s*(.+?)\s*=\s*(\S+?);?:?\s*$/.exec(line); return m ? [m[1], num(m[2])] : null; };
const xsOf = (line) => [...line.matchAll(/x(\d+) = (\S+?);/g)].map((m) => num(m[2]));

function parseIterations(lines) {
  // "Iteration n:" (or the labs' "Interation 0:") followed by "key = value" lines
  const its = [];
  for (const l of lines) {
    const h = /^\s*In?teration (\d+):/.exec(l);
    if (h) { its.push({ n: +h[1] }); continue; }
    const p = its.length && kv(l);
    if (p) its[its.length - 1][p[0]] = p[1];
  }
  return its;
}

function section(lines, start, ends) {
  const i = lines.findIndex((l) => l.includes(start));
  if (i < 0) return [];
  const out = [];
  for (const l of lines.slice(i + 1)) { if (ends.some((e) => l.includes(e))) break; out.push(l); }
  return out;
}

// --------------------------------------------------------------- the labs --
const iterCol = (its) => its.map((it) => it.n);

const LAB2_M = [
  [4.855, 1.239, 0.272, 0.258, 1.192],
  [1.491, 4.954, 0.124, 0.236, 0.256],
  [0.456, 0.285, 4.354, 0.254, 0.852],
  [0.412, 0.335, 0.158, 2.874, 0.862],
];
const LAB5_X = [0.092, 0.772, 1.385, 2.108, 2.938];
const LAB5_Y = [3.161, 1.357, -0.158, -0.129, -4.438];
const LAB6_X = [0.235, 0.569, 0.903, 1.237, 1.571, 1.906, 2.240, 2.574, 2.908];
const LAB6_Y = [1.082, 1.563, 2.182, 3.078, 3.816, 4.500, 5.443, 6.209, 7.082];
const F3 = (x) => 1 - 0.5 * x ** 2 * Math.log(x) + 0.3 * Math.sqrt(x);
const F4 = ([a, b]) => [Math.sin(a + 1.5) - b + 2.9, Math.cos(b - 2) + a];
const F7 = (x) => Math.sqrt(x + 1) * Math.log10(x + 1);
function newtonRoot() {
  let x = [0, 4];
  for (let k = 0; k < 40; k++) {
    const [f0, f1] = F4(x), J = [Math.cos(x[0] + 1.5), -1, 1, -Math.sin(x[1] - 2)], det = J[0] * J[3] - J[1] * J[2];
    x = [x[0] - (J[3] * f0 - J[1] * f1) / det, x[1] - (-J[2] * f0 + J[0] * f1) / det];
  }
  return x;
}

const LABS = [
  {
    id: 'lab2', tab: '2', file: '2-lab.cpp', name: 'Seidel iteration',
    intro: 'Solves a 4×4 linear system by iteration: each equation is divided by its diagonal element and rewritten as x<sub>i</sub> = Σ α<sub>ij</sub>x<sub>j</sub> + β<sub>i</sub>. The loop stops once the largest change between iterations falls below (1 − ‖α‖)/‖α‖ · ε, with ε = 0.001.',
    globals: true,
    bugs: [
      'Each new x<sub>i</sub> is computed from x<sub>i</sub> alone. The inner loop multiplies every coefficient α<sub>ij</sub> by <code>X_prev[i]</code> instead of the other unknowns, so this is neither Seidel nor Jacobi iteration. It settles on a point that isn\'t the solution (compare the JS reference below).',
      '<code>B[i] +=</code> sums into an uninitialised local array. WebAssembly\'s stack starts zeroed, as the native unoptimised build\'s does in practice, so the sums come out right. Compiled with <code>-O2</code>, the program loops forever.',
    ],
    hint: 'The page writes the table into the program\'s global <code>float M[4][5]</code> before <code>main()</code> starts. Systems whose iteration diverges are cut off after 6000 printed lines.',
    run: (inp) => runWasm('2', { M: inp.M.flat() }),
    parse(lines) {
      const i = lines.indexOf('Equations set:');
      const alpha = lines.slice(i + 1, i + 5).map((l) => l.trim().split(/\s+/).map(num));
      const its = [];
      let mode = '';
      for (const l of lines) {
        let m;
        if ((m = /^Iteration (\d+):/.exec(l))) { its.push({ n: +m[1], prev: [], cur: [], d: [] }); continue; }
        const it = its[its.length - 1];
        if (!it) continue;
        if (/Previous X/.test(l)) mode = 'prev';
        else if (/Current X/.test(l)) mode = 'cur';
        else if ((m = /^\s+x1 = /.exec(l))) it[mode] = xsOf(l);
        else if ((m = /Delta x(\d): (\S+)/.exec(l))) it.d[+m[1] - 1] = num(m[2]);
        else if ((m = /Current delta \(max\): (\S+);/.exec(l))) it.max = num(m[1]);
        else if ((m = /Pivotal misstep: (\S+)/.exec(l))) it.stop = num(m[1]);
      }
      const fin = lines.find((l) => l.includes('[FINAL RESULT]'));
      return { alpha, its, final: fin ? xsOf(fin) : null };
    },
    render(d, inp) {
      const A = inp.M.map((r) => r.slice(0, 4)), b = inp.M.map((r) => r[4]);
      const ref = gauss(A, b);
      const fin = d.final || d.its[d.its.length - 1]?.cur || [];
      const resid = fin.length ? Math.max(...A.map((r, i) => Math.abs(r.reduce((s, a, j) => s + a * fin[j], 0) - b[i]))) : NaN;
      const err = (x) => (ref ? Math.max(...x.map((v, j) => Math.abs(v - ref[j]))) : NaN);
      const n = d.its.length;
      const xsIt = [...iterCol(d.its), ...(d.final ? [n + 1] : [])];
      const comp = (j) => [...d.its.map((it) => it.cur[j]), ...(d.final ? [d.final[j]] : [])];
      const colors = ['--c1', '--c2', '--c3', '--c4'];
      return `<div class="stats">${stat(d.final ? n : `${n}+`, 'iterations')}${fin.map((v, j) => stat(fmt(v, 6), `x${j + 1}, lab result`)).join('')}${stat(fmt(d.its[0]?.stop, 6), 'stop threshold')}${stat(fmt(resid, 4), 'max |Ax − b| of the result', resid > 1e-3 ? 'bad' : 'good')}</div>
        ${ref ? `<div class="stats">${ref.map((v, j) => stat(fmt(v, 6), `x${j + 1}* ${JS} Gauss`)).join('')}${stat(fmt(err(fin), 4), `max |x − x*| ${JS}`, err(fin) > 1e-3 ? 'bad' : 'good')}</div>` : ''}
        <div class="charts two">
          ${chart({ title: 'Unknowns per iteration', x: 'iteration', y: 'x', intX: true, series: range(4, (j) => ({ name: `x${j + 1}`, xs: xsIt, ys: comp(j), color: colors[j] })), note: d.final ? `The last point (${n + 1}) is the final result.` : '' })}
          ${chart({ title: 'Convergence', x: 'iteration', y: 'change', intX: true, log: true, logToggle: true, series: [
            { name: 'max |Δx|', xs: iterCol(d.its), ys: d.its.map((it) => it.max), color: '--c1' },
            { name: 'stop threshold', xs: iterCol(d.its), ys: d.its.map((it) => it.stop), color: '--c5', dash: true },
            ...(ref ? [{ name: 'max |x − x*|', ref: true, xs: xsIt, ys: xsIt.map((_, k) => err(range(4, (j) => comp(j)[k]))), color: '--c7' }] : []),
          ] })}
        </div>
        ${table('Iterations', ['#', ...range(4, (j) => `x${j + 1}`), ...range(4, (j) => `|Δx${j + 1}|`), 'max |Δx|'], d.its.map((it) => [it.n, ...it.cur, ...it.d, it.max]), { note: 'x at the start of each iteration and the change it made' })}
        ${table('Coefficients α and β', ['row', 'α', 'α', 'α', 'β'], d.alpha.map((r, i) => [i + 1, ...r]), { note: 'as printed: the three off-diagonal α of each row, then β' })}`;
    },
  },
  {
    id: 'lab3', tab: '3', file: '3-lab.cpp', name: 'Chord method',
    intro: 'Finds the root of f(x) = 1 − 0.5·x²·ln x + 0.3·√x on [2, 3] with the chord (secant) method: the next x is where the chord through the last two points crosses zero. It stops when |x − x<sub>prev</sub>| ≤ ε = 10⁻⁶.',
    consts: ['A = 2', 'B = 3', 'eps = 1e-6'],
    run: () => runWasm('3'),
    parse(lines) {
      const its = parseIterations(lines);
      const fin = /\[FINAL RESULT\]: (\S+)/.exec(lines.join('\n'));
      return { its, root: fin ? num(fin[1]) : NaN };
    },
    render(d) {
      const its = d.its, last = its[its.length - 1];
      // the chord of step k passes through s[k] and s[k+1] and crosses zero at s[k+2]; s = [A, B, x1, x2, …]
      const s = [its[0].x, 3, ...its.slice(1).map((it) => it.x)];
      const fv = (k) => (k === 1 ? F3(3) : its[k === 0 ? 0 : k - 1]['f(x)']);
      const segments = range(Math.max(0, s.length - 2), (k) => ({ pts: [[s[k], fv(k)], [s[k + 1], fv(k + 1)], [s[k + 2], 0]], color: '--c2', dash: true }));
      const xs = linspace(1.95, 3.05, 200);
      return `<div class="stats">${stat(fmt(d.root, 7), 'root x')}${stat(fmt(last['f(x)'], 8), 'f(x) at the root')}${stat(its.length - 1, 'iterations')}${stat('10⁻⁶', 'ε')}</div>
        <div class="charts">${chart({ title: 'f(x) and the chords', x: 'x', y: 'f(x)', height: 280, series: [
          { name: 'f(x)', ref: true, xs, ys: xs.map(F3), color: '--c1' },
          { name: 'iterates x_k', xs: its.map((it) => it.x), ys: its.map((it) => it['f(x)']), color: '--c5', line: false, r: 3.5 },
        ], segments, note: 'Dashed lines are the chords the program steps along. Drag across the root near x ≈ 2.017 to zoom in.' })}</div>
        <div class="charts two" style="margin-top:16px">
          ${chart({ title: 'x per iteration', x: 'iteration', y: 'x', intX: true, series: [{ name: 'x', xs: iterCol(its), ys: its.map((it) => it.x), color: '--c1' }] })}
          ${chart({ title: 'Convergence', x: 'iteration', y: 'value', intX: true, log: true, logToggle: true, series: [
            { name: '|x − x_prev|', xs: iterCol(its).slice(1), ys: its.slice(1).map((it) => it['|x - x_prev|']), color: '--c3' },
            { name: '|f(x)|', xs: iterCol(its), ys: its.map((it) => Math.abs(it['f(x)'])), color: '--c4' },
          ] })}
        </div>
        ${table('Iterations', ['#', 'x', 'f(x)', '|x − x_prev|'], its.map((it) => [it.n, it.x, it['f(x)'], it['|x - x_prev|']]))}`;
    },
  },
  {
    id: 'lab41', tab: '4.1', file: '4.1-lab.cpp', wasm: '4.1', name: 'Newton method',
    intro: 'Solves the nonlinear system sin(x₁ + 1.5) − x₂ + 2.9 = 0, cos(x₂ − 2) + x₁ = 0 from (0, 4) with Newton steps x ← x − J⁻¹F(x), until the largest change is below ε = 10⁻⁶ or after 21 iterations.',
    consts: ['X1 = 0', 'X2 = 4', 'eps = 1e-6', 'ITERATIONS = 20'],
    bugs: [
      'The Jacobian is filled in transposed (<code>J[1] = 1</code>, <code>J[2] = -1</code>; ∂F₁/∂x₂ is −1 and ∂F₂/∂x₁ is 1). <code>Jackobian_min1</code> then swaps the off-diagonal entries and divides by the determinant, which isn\'t the inverse ([d −b; −c a]/det). The steps aren\'t Newton steps, so the iterates run away until the iteration cap. The method is also not the "simplified" variant, because it recomputes the Jacobian every step.',
    ],
    run: () => runWasm('4.1'),
    parse: (lines) => ({ its: parseIterations(lines) }),
    render: (d) => systemView(d),
  },
  {
    id: 'lab42', tab: '4.2', file: '4.2-lab.cpp', wasm: '4.2', name: 'Steepest descent',
    intro: 'Minimises Φ(x) = F₁² + F₂² for the same system as lab 4.1, from (0, 4) with step a = 0.3, until the largest change is below ε = 10⁻⁶ or after 21 iterations.',
    consts: ['X1 = 0', 'X2 = 4', 'eps = 1e-6', 'a = 0.3', 'ITERATIONS = 20'],
    bugs: [
      'The loop replaces x with the gradient expressions themselves (<code>X[0] = ∂Φ/∂x₁</code>, then <code>X[1] = ∂Φ/∂x₂</code> using the already overwritten <code>X[0]</code>) instead of stepping x − a·∇Φ(x). The step <code>a</code> is declared but never used, and the expressions don\'t match the derivatives of Φ (for example −5.8·cos(x₁ + 1.5), where the derivative has +5.8). So the iterates grow until the iteration cap.',
    ],
    run: () => runWasm('4.2'),
    parse: (lines) => ({ its: parseIterations(lines) }),
    render: (d) => systemView(d),
  },
  {
    id: 'lab5', tab: '5', file: '5-lab.cpp', name: 'Interpolation',
    intro: 'Interpolates a five-point table with the Lagrange polynomial L(x) and the Newton polynomial N(x) built from divided differences, prints both polynomials and evaluates them at x = x₁ + x₂ (the second and third nodes).',
    globals: true,
    bugs: [
      'The "finite differences" step prints the shifted nodes x<sub>i</sub> + k·h with h = 1/SIZE = 0.2, not differences of y.',
      '<code>float res</code> in <code>LxFunction</code> isn\'t initialised before <code>res += …</code>. It happens to start at 0 in WebAssembly and in the native unoptimised build.',
    ],
    hint: 'The page writes the table into the program\'s global <code>float X[5]</code> and <code>Y[5]</code> before <code>main()</code> starts. The curves are drawn from the polynomial coefficients the program prints.',
    run: (inp) => runWasm('5', { X: inp.X, Y: inp.Y }),
    parse(lines) {
      const terms = (block) => [...block.join(' ').matchAll(/(-?\d+\.\d+|-?nan|-?inf)((?:\(x - -?[\d.]+\))*)/g)]
        .map((m) => ({ c: num(m[1]), roots: [...m[2].matchAll(/\(x - (-?[\d.]+)\)/g)].map((r) => +r[1]) }));
      const ends = ['P(x):', 'L(x):', 'L(x1', 'Finite', 'Divided', 'N(x):', 'N(x1'];
      const P = section(lines, 'P(x):', ends.filter((e) => e !== 'P(x):')).filter((l) => l.trim()).map((l) => terms([l.replace(/^x\d:\s*/, '')])[0]);
      const L = terms(section(lines, 'L(x):', ['L(x1']));
      const N = terms(section(lines, 'N(x):', ['N(x1']));
      const tri = (name) => section(lines, name, ends.filter((e) => e !== name)).filter((l) => l.trim()).map((l) => l.trim().split(/\s+/).map(num));
      const val = (name) => { const m = new RegExp(`${name.replace(/[()+]/g, '\\$&')} = (\\S+)`).exec(lines.join('\n')); return m ? num(m[1]) : NaN; };
      return { P, L, N, fin: tri('Finite differences:'), div: tri('Divided differences:'), Lv: val('L(x1 + x2)'), Nv: val('N(x1 + x2)') };
    },
    render(d, inp) {
      const ev = (ts) => (x) => ts.reduce((s, t) => s + t.c * t.roots.reduce((p, r) => p * (x - r), 1), 0);
      const lo = Math.min(...inp.X), hi = Math.max(...inp.X), pad = (hi - lo) * 0.06;
      const xs = linspace(lo - pad, hi + pad, 240), xq = inp.X[1] + inp.X[2];
      const order = inp.X.map((x, i) => i).sort((a, b) => inp.X[a] - inp.X[b]);
      return `<div class="stats">${stat(fmt(xq, 4), 'x = x₁ + x₂')}${stat(fmt(d.Lv, 4), 'L(x₁ + x₂)')}${stat(fmt(d.Nv, 4), 'N(x₁ + x₂)')}${stat(fmt(Math.abs(d.Lv - d.Nv), 4), '|L − N|', Math.abs(d.Lv - d.Nv) < 1e-3 ? 'good' : 'bad')}</div>
        <div class="charts">${chart({ title: 'Interpolating polynomials and the table', x: 'x', y: 'y', height: 290, series: [
          { name: 'L(x)', xs, ys: xs.map(ev(d.L)), color: '--c1', width: 2.6 },
          { name: 'N(x)', xs, ys: xs.map(ev(d.N)), color: '--c2', dash: true },
          { name: 'table points', xs: order.map((i) => inp.X[i]), ys: order.map((i) => inp.Y[i]), color: '--c5', line: false, r: 4 },
        ], marks: [{ x: xq, y: d.Lv, label: `x₁ + x₂ → ${fmt(d.Lv, 4)}`, color: '--c3' }], note: 'L(x) and N(x) are evaluated from the coefficients the program prints (6 decimals), so they agree to about that precision.' })}</div>
        <div class="section-h"><h2>Newton polynomial <small>· N(x) as printed</small></h2></div><div class="formula">${esc(d.N.map((t) => `${fmt(t.c, 6)}${t.roots.map((r) => `(x − ${r})`).join('')}`).join(' + '))}</div>
        <div class="two-col" style="margin-top:4px">
          ${table('Divided differences', ['order', ...range(5, (i) => `#${i}`)], d.div.map((r, k) => [k, ...range(k, () => ''), ...r]))}
          ${table('"Finite differences"', ['row', ...range(5, (i) => `k = ${i}`)], d.fin.map((r, k) => [k, ...r]), { note: 'nodes + k·h, as printed' })}
        </div>
        ${table('Lagrange basis', ['i', 'x_i', 'y_i', '1 / Π(x_i − x_j)', 'y_i / Π(x_i − x_j)'], d.P.map((t, i) => [i, inp.X[i], inp.Y[i], t.c, d.L[i]?.c]))}`;
    },
  },
  {
    id: 'lab6', tab: '6', file: '6-lab.cpp', name: 'Least squares',
    intro: 'Approximates a nine-point table with two model functions, y = a/x² + b/x + c and y = a/x + b·eˣ + c. It prints the normal system of the least-squares method for each, then the sum of squared deviations S and δ = √S for the coefficients set in <code>main()</code>.',
    consts: ['SIZE = 9', 'a1 = 0.1177, b1 = 1.7591, c1 = 1.681', 'a2 = −0.4091, b2 = −0.5904, c2 = 12.2568'],
    bugs: [
      '<code>main()</code> plugs in fixed coefficients instead of solving the normal systems it prints, so both fits are poor (S = 170 and 348). The dashed JS curves are the actual least-squares fits.',
      '<code>func2Sys</code> adds <code>Y[i]</code> to <code>sumY</code> twice and adds both y·x and y/x to <code>sumYX</code>, so the printed system for function 2 is wrong too.',
    ],
    run: () => runWasm('6'),
    parse(lines) {
      const sys = (name) => section(lines, `${name} (normal system):`, ['Function']).filter((l) => l.includes('=')).map((l) => [...l.matchAll(/(-?\d+\.\d+)/g)].map((m) => +m[1]));
      const fn = (name) => {
        const i = lines.findIndex((l) => l.trim() === `${name}:`);
        const co = [...(lines[i + 1] || '').matchAll(/(-?\d+\.\d+)/g)].map((m) => +m[1]);
        return { a: co[0], b: co[1], c: co[2], S: kv(lines[i + 2] || '')?.[1], del: kv(lines[i + 3] || '')?.[1] };
      };
      return { sys1: sys('Function 1'), sys2: sys('Function 2'), f1: fn('Function 1'), f2: fn('Function 2') };
    },
    render(d) {
      const B1 = [(x) => x ** -2, (x) => 1 / x, () => 1], B2 = [(x) => 1 / x, Math.exp, () => 1];
      const ls1 = leastSquares(LAB6_X, LAB6_Y, B1), ls2 = leastSquares(LAB6_X, LAB6_Y, B2);
      const f = (B, c) => (x) => B.reduce((s, g, j) => s + c[j] * g(x), 0);
      const S = (fx) => LAB6_X.reduce((s, x, i) => s + (fx(x) - LAB6_Y[i]) ** 2, 0);
      const solved1 = gauss(d.sys1.map((r) => r.slice(0, 3)), d.sys1.map((r) => r[3]));
      const xs = linspace(0.2, 2.95, 220);
      const lab1 = f(B1, [d.f1.a, d.f1.b, d.f1.c]), lab2 = f(B2, [d.f2.a, d.f2.b, d.f2.c]);
      const card = (name, form, fl, ls, Sls, sysRows, solved) => `<div class="panel" style="margin:0;padding:14px"><h2 style="margin-bottom:8px">${name} <small>· ${form}</small></h2>
        <div class="stats" style="margin-bottom:10px">${stat(fmt(fl.S, 4), 'S, lab coefficients', 'bad')}${stat(fmt(fl.del, 4), 'δ = √S')}${stat(fmt(Sls, 4), `S, least squares ${JS}`, 'good')}</div>
        <div class="table-wrap"><table><thead><tr><th class="num"></th><th class="num">a</th><th class="num">b</th><th class="num">c</th></tr></thead><tbody>
          <tr><td>lab (main)</td>${[fl.a, fl.b, fl.c].map((v) => `<td class="num">${fmt(v, 4)}</td>`).join('')}</tr>
          ${solved ? `<tr><td>printed system solved ${JS}</td>${solved.map((v) => `<td class="num">${fmt(v, 4)}</td>`).join('')}</tr>` : ''}
          <tr><td>least squares ${JS}</td>${ls.map((v) => `<td class="num">${fmt(v, 4)}</td>`).join('')}</tr></tbody></table></div>
        ${table(`${name}: normal system`, ['row', '· a', '· b', '· c', '= rhs'], sysRows.map((r, i) => [i + 1, ...r]), { note: 'as printed' })}</div>`;
      return `<div class="charts">${chart({ title: 'Model functions against the table', x: 'x', y: 'y', height: 300, view: { y: [-2, 12] }, series: [
          { name: 'table points', xs: LAB6_X, ys: LAB6_Y, color: '--c5', line: false, r: 4 },
          { name: 'function 1 (lab)', xs, ys: xs.map(lab1), color: '--c1' },
          { name: 'function 2 (lab)', xs, ys: xs.map(lab2), color: '--c2' },
          { name: 'function 1, LS', ref: true, xs, ys: xs.map(f(B1, ls1)), color: '--c1', dash: true },
          { name: 'function 2, LS', ref: true, xs, ys: xs.map(f(B2, ls2)), color: '--c2', dash: true },
        ], note: 'The table is the <code>X</code>, <code>Y</code> arrays from the source. Function 1 shoots up near x = 0, so the default view is clipped; use Full range to see all of it.' })}</div>
        <div class="algo-cards" style="margin-top:16px">
          ${card('Function 1', 'a/x² + b/x + c', d.f1, ls1, S(f(B1, ls1)), d.sys1, solved1)}
          ${card('Function 2', 'a/x + b·eˣ + c', d.f2, ls2, S(f(B2, ls2)), d.sys2, null)}
        </div>`;
    },
  },
  {
    id: 'lab7', tab: '7', file: '7-lab.cpp', name: 'Integration',
    intro: 'Integrates f(x) = √(x + 1)·lg(x + 1) over [0.1, 1.1] with the trapezoid rule and Simpson\'s rule, each with step h = 0.125 (8 steps) and 2h, and estimates the error by Runge\'s rule: |I(h) − I(2h)|/3 for trapezoids and /15 for Simpson.',
    consts: ['A = 0.1', 'B = 1.1', 'h = 0.125', 'n = 8'],
    bugs: [
      'Simpson\'s rule multiplies by <code>h * 0.6667</code> instead of 2h/3. That shifts I(h) by about 1.3·10⁻⁵, far more than the 10⁻⁷ error estimate it prints (compare the JS integral).',
    ],
    run: () => runWasm('7'),
    parse(lines) {
      const block = (name) => Object.fromEntries(section(lines, name, ['method:']).map(kv).filter(Boolean));
      return { trap: block('Trapeze method:'), simp: block("Simpson's method:") };
    },
    render(d) {
      const exact = simpson(F7, 0.1, 1.1, 20000);
      const xs = linspace(0.05, 1.15, 160), nodes = linspace(0.1, 1.1, 9);
      const r = [['Trapezoid', d.trap], ['Simpson', d.simp]];
      const hs = [0.125, 0.25];
      return `<div class="stats">${stat(fmt(d.trap['I(h)'], 6), 'trapezoid I(h)')}${stat(fmt(d.trap['I-I^tp(h)'], 6), 'trapezoid error estimate')}${stat(fmt(d.simp['I(h)'], 8), 'Simpson I(h)')}${stat(d.simp['I-I^tp(h)'].toExponential(2), 'Simpson error estimate')}${stat(fmt(exact, 8), `integral ${JS}`)}</div>
        <div class="charts two">
          ${chart({ title: 'f(x) with the trapezoid strips, h = 0.125', x: 'x', y: 'f(x)', height: 260, series: [{ name: 'f(x)', ref: true, xs, ys: xs.map(F7), color: '--c1' }], areas: [{ xs: nodes, ys: nodes.map(F7), color: '--c7' }], note: 'f is evaluated by this page to draw the curve; the integrals come from the program.' })}
          ${chart({ title: 'Integral against the step size', x: 'step h', y: 'I', height: 260, xPad: 0.15, dots: true, series: [
            { name: 'trapezoid', xs: hs, ys: [d.trap['I(h)'], d.trap['I(2h)']], color: '--c2' },
            { name: 'Simpson', xs: hs, ys: [d.simp['I(h)'], d.simp['I(2h)']], color: '--c3' },
            { name: 'integral', ref: true, xs: hs, ys: [exact, exact], color: '--c5', dash: true },
          ] })}
        </div>
        ${table('Results', ['method', 'I(h), h = 0.125', 'I(2h), 2h = 0.25', 'error estimate', `|I(h) − integral| ${JS}`], r.map(([n, v]) => [n, v['I(h)'], v['I(2h)'], v['I-I^tp(h)'], Math.abs(v['I(h)'] - exact)]))}`;
    },
  },
  {
    id: 'lab8', tab: '8', file: '8-lab.cpp', name: 'Runge–Kutta and Euler',
    intro: 'Solves y′ = 0.5·(x + 1)·eˣ·y² − x·y, y(0) = 1 on [0, 2]. It first checks the step h = 0.1 by comparing one 2h step with two h steps (doubling h while they agree within ε = 10⁻⁴), then runs fourth-order Runge–Kutta and Euler, and compares RK with the exact solution.',
    consts: ['A = 0, B = 2', 'h = 0.1', 'eps = 1e-4', 'x0 = 0, y0 = 1'],
    bugs: [
      'The fourth Runge–Kutta stage uses x + h/2 instead of x + h.',
      'Euler adds h·(the whole RK step result) instead of h·f(x, y).',
      'The exact comparison pairs y(x<sub>i</sub>) with the RK value at x<sub>i</sub> + h.',
      '<code>RK_Yvalues[20]</code> is written at index 20 on the last step (the float x loop runs 20 times), one past the end of the array.',
    ],
    hint: 'The exact solution has a pole between x = 0.9 and 1.0, so every method blows up there and the output turns to inf and nan.',
    run: () => runWasm('8'),
    parse(lines) {
      const blocks = (name) => parseIterations(section(lines, name, ['~~~~', 'Max delta']));
      const max = (name) => { const s = section(lines, name, ['~~~~']).find((l) => l.startsWith('Max delta')); return s ? num(s.split('=')[1]) : NaN; };
      const search = [];
      for (const l of section(lines, '[STEP SEARCHING]', ['[AFTER STEP'])) {
        const p = kv(l);
        if (!p) continue;
        if (p[0] === 'y1') search.push({});
        search[search.length - 1][p[0]] = p[1];
      }
      const after = Object.fromEntries(section(lines, '[AFTER STEP SEARCHING]', ['~~~~']).map(kv).filter(Boolean));
      return { search, after, rk: blocks('[RUNGE-KUTTA(IV) METHOD]'), eu: blocks('[EULER METHOD]'), ex: blocks('[EXACT CALCULATING]'), maxRk: max('[RUNGE-KUTTA(IV) METHOD]'), maxEu: max('[EULER METHOD]'), maxEx: max('[EXACT CALCULATING]') };
    },
    render(d) {
      const h = d.after.h0;
      const step = (b) => b.map((it) => it.xi + h);
      const rows = (b, k2, k3) => b.map((it) => [it.n, it.xi, it.yi, it[k2], it[k3]]);
      return `<div class="stats">${stat(fmt(h, 4), 'step h₀')}${stat(d.after.n, 'steps n')}${stat(fmt(d.after['|y2 - ~y2|'], 6), '|y₂ − ỹ₂| at h₀')}${stat(fmt(d.maxRk), 'RK max |y − ỹ|', 'bad')}${stat(fmt(d.maxEu), 'Euler max |y − ỹ|', 'bad')}${stat(fmt(d.maxEx), 'max |y − y_RK|', 'bad')}</div>
        <div class="charts">${chart({ title: 'Solutions', x: 'x', y: 'y', height: 300, view: { y: [-8, 10] }, series: [
          { name: 'exact', xs: d.ex.map((it) => it.xi), ys: d.ex.map((it) => it.yi), color: '--c3', width: 2.4 },
          { name: 'Runge–Kutta', xs: step(d.rk), ys: d.rk.map((it) => it.yi), color: '--c1' },
          { name: 'Euler', xs: step(d.eu), ys: d.eu.map((it) => it.yi), color: '--c2' },
        ], note: 'RK and Euler values are plotted at x<sub>i</sub> + h, the point each step reaches. The y axis is clipped to −8…10 by default; Full range shows the blow-up.' })}</div>
        <div class="charts two" style="margin-top:16px">
          ${chart({ title: 'Step check: |y(2h) − ỹ(2h)|', x: 'x', y: 'difference', log: true, logToggle: true, series: [
            { name: 'Runge–Kutta', xs: d.rk.map((it) => it.xi), ys: d.rk.map((it) => it['|yi - ~yi|']), color: '--c1' },
            { name: 'Euler', xs: d.eu.map((it) => it.xi), ys: d.eu.map((it) => it['|yi - ~yi|']), color: '--c2' },
          ] })}
          ${chart({ title: 'Error against the exact solution', x: 'x', y: '|y − y_RK|', log: true, logToggle: true, series: [{ name: '|y − y_RK|', xs: d.ex.map((it) => it.xi), ys: d.ex.map((it) => it['|yi - ~yi|']), color: '--c4' }] })}
        </div>
        <div data-copyable data-title="Iterations" style="margin-top:4px"><div class="section-h"><h2>Iterations</h2>${actions('table')}</div>
          <div class="tabs" data-tabs>${['Runge–Kutta', 'Euler', 'Exact', 'Step search'].map((t, i) => `<button type="button" class="tab" data-tab="${i}" aria-selected="${i === 0}">${t}</button>`).join('')}</div>
          ${[
            [['#', 'x_i', 'y_i', '~y_i', '|y_i − ~y_i|'], rows(d.rk, '~yi', '|yi - ~yi|')],
            [['#', 'x_i', 'y_i', '~y_i', '|y_i − ~y_i|'], rows(d.eu, '~yi', '|yi - ~yi|')],
            [['#', 'x_i', 'y_i exact', 'y_i RK', '|y_i − y_RK|'], rows(d.ex, 'yi_RK', '|yi - ~yi|')],
            [['pass', 'y₁', 'y₂', '~y₂', '|y₂ − ~y₂|', 'h₀'], [...d.search.map((s, i) => [i + 1, s.y1, s.y2, s['~y2'], s['|y2 - ~y2|'], s.h0]), ['after', d.after.y1, d.after.y2, d.after['~y2'], d.after['|y2 - ~y2|'], h]]],
          ].map(([head, body], i) => `<div class="table-wrap" data-pane="${i}"${i ? ' hidden' : ''}><table><thead><tr>${head.map((x) => `<th class="num">${x}</th>`).join('')}</tr></thead><tbody>${body.map((r) => `<tr>${r.map((c) => `<td class="num">${typeof c === 'number' ? fmt(c, 6) : c}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`).join('')}
        </div>`;
    },
  },
];

function systemView(d) {
  const its = d.its, last = its[its.length - 1], root = newtonRoot();
  const converged = last.delta < 1e-6;
  return `<div class="stats">${stat(its.length - 1, 'iterations')}${stat(fmt(last.x1, 4), 'last x₁')}${stat(fmt(last.x2, 4), 'last x₂')}${stat(fmt(last.delta, 4), 'last max |Δx|', converged ? 'good' : 'bad')}${stat(converged ? 'converged' : 'iteration cap', 'stopped by', converged ? 'good' : 'bad')}${stat(`(${fmt(root[0], 5)}, ${fmt(root[1], 5)})`, `root ${JS} Newton`)}</div>
    <div class="charts two">
      ${chart({ title: 'Convergence: max |Δx| per iteration', x: 'iteration', y: 'max |Δx|', intX: true, log: true, logToggle: true, series: [
        { name: 'max |Δx|', xs: iterCol(its).slice(1), ys: its.slice(1).map((it) => it.delta), color: '--c1' },
        { name: 'ε', xs: iterCol(its).slice(1), ys: its.slice(1).map(() => 1e-6), color: '--c5', dash: true },
      ] })}
      ${chart({ title: 'Residual |F(x)| per iteration', x: 'iteration', y: 'max |F(x)|', intX: true, log: true, logToggle: true, series: [{ name: 'max |F(x)|', ref: true, xs: iterCol(its), ys: its.map((it) => Math.max(...F4([it.x1, it.x2]).map(Math.abs))), color: '--c4' }], note: 'F is evaluated by this page at each printed iterate.' })}
    </div>
    <div class="charts" style="margin-top:16px">${chart({ title: 'Iterates x₁, x₂', x: 'iteration', y: 'x', intX: true, series: [
      { name: 'x₁', xs: iterCol(its), ys: its.map((it) => it.x1), color: '--c2' },
      { name: 'x₂', xs: iterCol(its), ys: its.map((it) => it.x2), color: '--c3' },
      { name: 'x₁*', ref: true, xs: iterCol(its), ys: its.map(() => root[0]), color: '--c2', dash: true, dots: false, width: 1.2 },
      { name: 'x₂*', ref: true, xs: iterCol(its), ys: its.map(() => root[1]), color: '--c3', dash: true, dots: false, width: 1.2 },
    ], note: 'Drag over the first iterations to zoom in; the dashed lines are the system\'s root.' })}</div>
    ${table('Iterations', ['#', 'x₁', 'x₂', 'max |Δx|'], its.map((it) => [it.n, it.x1, it.x2, it.delta]))}`;
}

// ------------------------------------------------------------------ copying --
const ACTION_LABEL = { png: 'Copy image', csv: 'Copy CSV', table: 'Copy table', results: 'Copy results', output: 'Copy output' };
const ACTION_TIP = { png: 'Copy as a PNG image (downloads it if the clipboard is unavailable)', csv: 'Copy the plotted data as CSV', table: 'Copy the table as tab-separated text, pastes into spreadsheets', results: 'Copy the inputs and results as text', output: 'Copy the program\'s console output as text' };
const actionButtons = (...kinds) => kinds.map((k) => `<button type="button" class="mini" data-act="${k}" title="${ACTION_TIP[k]}">${ACTION_LABEL[k]}</button>`).join('');
const actions = (...kinds) => `<span class="copy-actions">${actionButtons(...kinds)}</span>`;

// Render an on-screen SVG to PNG: computed colours and fonts are inlined so the image
// matches the current theme, with the caption drawn above it.
function svgPng(svg, title, legend = [], scale = 2) {
  const box = svg.viewBox.baseVal;
  const w = (box && box.width) || +svg.getAttribute('width') || svg.clientWidth;
  const h = (box && box.height) || +svg.getAttribute('height') || svg.clientHeight;
  const clone = svg.cloneNode(true);
  const src = [svg, ...svg.querySelectorAll('*')], dst = [clone, ...clone.querySelectorAll('*')];
  src.forEach((el, i) => {
    const cs = getComputedStyle(el), c = dst[i];
    for (const a of ['fill', 'stroke', 'stroke-opacity', 'fill-opacity', 'stroke-width', 'opacity', 'visibility']) c.setAttribute(a, cs.getPropertyValue(a));
    if (el.tagName === 'text') for (const a of ['font-family', 'font-size', 'font-weight']) c.setAttribute(a, cs.getPropertyValue(a));
    c.removeAttribute('style'); c.removeAttribute('class');
  });
  clone.querySelectorAll('title, .marks').forEach((t) => t.remove());
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width', w); clone.setAttribute('height', h);
  if (!clone.getAttribute('viewBox')) clone.setAttribute('viewBox', `0 0 ${w} ${h}`);
  const css = getComputedStyle(document.documentElement);
  const head = title ? 30 : 0, pad = 14;
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml' }));
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => {
      const c = document.createElement('canvas');
      c.width = (w + pad * 2) * scale; c.height = (h + head + pad * 2) * scale;
      const ctx = c.getContext('2d');
      ctx.scale(scale, scale);
      ctx.fillStyle = css.getPropertyValue('--panel').trim();
      ctx.fillRect(0, 0, w + pad * 2, h + head + pad * 2);
      const family = getComputedStyle(document.body).fontFamily;
      if (title) {
        ctx.fillStyle = css.getPropertyValue('--ink').trim();
        ctx.font = `600 14px ${family}`;
        ctx.fillText(title, pad, pad + 14);
        ctx.font = `12px ${family}`;
        let x = w + pad;
        for (const { name, color } of [...legend].reverse()) {
          x -= ctx.measureText(name).width;
          ctx.fillStyle = css.getPropertyValue('--muted').trim();
          ctx.fillText(name, x, pad + 13);
          x -= 19;
          ctx.fillStyle = css.getPropertyValue(color).trim();
          ctx.fillRect(x, pad + 7, 14, 3);
          x -= 14;
        }
      }
      ctx.drawImage(img, pad, pad + head, w, h);
      URL.revokeObjectURL(url);
      c.toBlob((b) => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('SVG render failed')); };
    img.src = url;
  });
}

const csvCell = (v) => (/[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v));
function chartCsv(spec) {
  const s0 = spec.series[0];
  const shared = spec.series.every((s) => s.xs.length === s0.xs.length && s.xs.every((x, i) => x === s0.xs[i]));
  if (shared) return [[spec.x, ...spec.series.map((s) => s.name || spec.y)], ...s0.xs.map((x, i) => [x, ...spec.series.map((s) => s.ys[i])])].map((r) => r.map(csvCell).join(',')).join('\n');
  return [['series', spec.x, spec.y], ...spec.series.flatMap((s) => s.xs.map((x, i) => [s.name, x, s.ys[i]]))].map((r) => r.map(csvCell).join(',')).join('\n');
}
const tableTsv = (table) => [...table.rows].map((r) => [...r.cells].map((c) => c.innerText.trim().replace(/\s+/g, ' ')).join('\t')).join('\n');

function resultsText(sec) {
  const lab = LABS.find((l) => `panel-${l.id}` === sec.id);
  const lines = [`Lab ${lab.tab}: ${lab.name} (${lab.file})`];
  const inp = state.get(lab.id).inputs;
  if (lab.id === 'lab2') lines.push('System [A | b]:', ...inp.M.map((r) => `  ${r.join('\t')}`));
  else if (lab.id === 'lab5') lines.push(`X: ${inp.X.join(', ')}`, `Y: ${inp.Y.join(', ')}`);
  else lines.push(`Constants: ${lab.consts.join(', ')}`);
  const out = $('.out', sec);
  out.querySelectorAll('.stats').forEach((st) => {
    const group = st.closest('.panel:not([role])')?.querySelector('h2')?.innerText.trim();
    if (group) lines.push(`${group}:`);
    st.querySelectorAll('.stat').forEach((x) => lines.push(`${group ? '  ' : ''}${$('span', x).innerText.trim()}: ${$('b', x).innerText.trim()}`));
  });
  return lines.join('\n');
}

async function copyAction(btn) {
  const kind = btn.dataset.act, box = btn.closest('[data-copyable]');
  const done = (text) => { btn.textContent = text; setTimeout(() => (btn.textContent = ACTION_LABEL[kind]), 1400); };
  const name = (box?.dataset.title || 'chart').toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
  try {
    if (kind === 'png') {
      const svg = box.querySelector('.plot svg');
      const st = charts.get(box.dataset.chart), spec = st?.spec;
      const legend = spec && spec.series.length > 1 ? spec.series.filter((x, i) => !st.hidden.has(i)).map((x) => ({ name: x.name, color: x.color })) : [];
      try {
        // pass the promise straight to ClipboardItem so Safari keeps the user gesture
        await navigator.clipboard.write([new ClipboardItem({ 'image/png': svgPng(svg, box.dataset.title, legend) })]);
        return done('Copied');
      } catch {
        const a = document.createElement('a');
        a.href = URL.createObjectURL(await svgPng(svg, box.dataset.title, legend));
        a.download = `${name}.png`;
        a.click();
        setTimeout(() => URL.revokeObjectURL(a.href), 1000);
        return done('Downloaded');
      }
    }
    const text = kind === 'csv' ? chartCsv(charts.get(box.dataset.chart).spec)
      : kind === 'table' ? tableTsv(box.querySelector('.table-wrap:not([hidden]) table'))
      : kind === 'output' ? state.get(btn.closest('section.panel').id.slice(6)).output
      : resultsText(btn.closest('section.panel'));
    await navigator.clipboard.writeText(text);
    done('Copied');
  } catch (e) { console.error(e); done('Failed'); }
}
document.addEventListener('click', (e) => { const b = e.target.closest('button[data-act]'); if (b) copyAction(b); });
document.addEventListener('click', (e) => {
  const t = e.target.closest('[data-tabs] .tab');
  if (!t) return;
  const box = t.closest('[data-copyable]');
  box.querySelectorAll('[data-tabs] .tab').forEach((x) => x.setAttribute('aria-selected', x === t));
  box.querySelectorAll('[data-pane]').forEach((p) => (p.hidden = p.dataset.pane !== t.dataset.tab));
});

// ------------------------------------------------------------------ inputs --
function inputsHtml(lab) {
  if (lab.id === 'lab2') {
    const head = ['', 'x₁', 'x₂', 'x₃', 'x₄', 'b'].map((h) => `<span class="h">${h}</span>`).join('');
    const cells = LAB2_M.map((r, i) => `<span class="h">${i + 1}</span>${r.map((v, j) => `<input name="m${i}${j}" inputmode="decimal" value="${v}" aria-label="row ${i + 1}, ${j < 4 ? `x${j + 1}` : 'b'}">`).join('')}`).join('');
    return `<div class="grid-scroll"><div class="grid-in" style="grid-template-columns:22px repeat(5, minmax(52px, 1fr));max-width:560px">${head}${cells}</div></div>`;
  }
  if (lab.id === 'lab5') {
    const row = (k, vals) => `<span class="h">${k}</span>${vals.map((v, i) => `<input name="${k}${i}" inputmode="decimal" value="${v}" aria-label="${k}${i}">`).join('')}`;
    return `<div class="grid-scroll"><div class="grid-in" style="grid-template-columns:22px repeat(5, minmax(52px, 1fr));max-width:520px">${['', ...range(5, (i) => `i = ${i}`)].map((h) => `<span class="h">${h}</span>`).join('')}${row('X', LAB5_X)}${row('Y', LAB5_Y)}</div></div>`;
  }
  return `<div class="form-row"><span class="hint" style="margin:0">Constants in the source:</span><ul class="consts">${lab.consts.map((c) => `<li>${esc(c)}</li>`).join('')}</ul></div>`;
}

function readInputs(lab, form) {
  const val = (name) => {
    const el = form.elements[name], v = Number(String(el.value).trim().replace(',', '.'));
    el.classList.toggle('bad', !Number.isFinite(v) || String(el.value).trim() === '');
    return v;
  };
  if (lab.id === 'lab2') {
    const M = range(4, (i) => range(5, (j) => val(`m${i}${j}`)));
    if (M.flat().some((v) => !Number.isFinite(v))) return { error: 'Every coefficient must be a number.' };
    if (M.some((r, i) => r[i] === 0)) return { error: 'The diagonal coefficients must be nonzero: the program divides each row by them.' };
    return { M };
  }
  if (lab.id === 'lab5') {
    const X = range(5, (i) => val(`X${i}`)), Y = range(5, (i) => val(`Y${i}`));
    if ([...X, ...Y].some((v) => !Number.isFinite(v))) return { error: 'Every table value must be a number.' };
    if (new Set(X.map((x) => Math.fround(x))).size < 5) return { error: 'The nodes X must be distinct: the program divides by their differences.' };
    return { X, Y };
  }
  return {};
}

// --------------------------------------------------------------- framework --
const nav = $('#nav'), main = $('#labs');
const state = new Map(); // lab id -> { ran, running, output, inputs }

function build() {
  for (const lab of LABS) {
    const b = document.createElement('button');
    b.type = 'button'; b.role = 'tab'; b.dataset.lab = lab.id;
    b.setAttribute('aria-controls', `panel-${lab.id}`);
    b.innerHTML = `<span class="lang cpp">C++</span>Lab ${lab.tab}<span class="long muted" style="font-weight:400">· ${esc(lab.name)}</span>`;
    b.onclick = () => select(lab.id, true);
    nav.appendChild(b);

    const sec = document.createElement('section');
    sec.className = 'panel'; sec.id = `panel-${lab.id}`; sec.role = 'tabpanel'; sec.hidden = true;
    sec.innerHTML = `<div class="lab-head"><span class="lang cpp">C++ → Wasm</span><h3>Lab ${lab.tab}: ${esc(lab.name)}</h3><a href="${SRC}/${lab.file}" target="_blank" rel="noopener">View source</a></div>
      <p class="intro">${lab.intro}</p>
      ${lab.bugs ? `<ul class="bugs">${lab.bugs.map((t) => `<li><b>Lab bug, kept as is:</b> ${t}</li>`).join('')}</ul>` : ''}
      <form class="params" novalidate>${inputsHtml(lab)}<div class="form-row"><button class="run" type="submit">Run</button>${lab.globals ? '<button class="ghost" type="button" data-reset-inputs>Reset to the lab\'s data</button>' : ''}</div></form>
      <p class="hint">${lab.globals ? lab.hint : `This lab reads no input: its data are compile-time constants, so Run re-runs the compiled program as is.${lab.hint ? ` ${lab.hint}` : ''}`}</p>
      <div class="out" style="margin-top:18px"></div>
      <details class="src"><summary>Source: <code>sripts/${lab.file}</code></summary><div data-copyable data-title="${lab.file}"><pre class="console"><code>Loading…</code></pre></div></details>`;
    main.appendChild(sec);
    state.set(lab.id, { ran: false, running: false, output: '', inputs: {} });

    const form = $('form', sec);
    form.addEventListener('submit', (e) => { e.preventDefault(); run(lab); });
    $('[data-reset-inputs]', sec)?.addEventListener('click', () => {
      const def = lab.id === 'lab2' ? Object.fromEntries(LAB2_M.flatMap((r, i) => r.map((v, j) => [`m${i}${j}`, v])))
        : Object.fromEntries([...LAB5_X.map((v, i) => [`X${i}`, v]), ...LAB5_Y.map((v, i) => [`Y${i}`, v])]);
      for (const [k, v] of Object.entries(def)) { form.elements[k].value = v; form.elements[k].classList.remove('bad'); }
      run(lab);
    });
    const det = $('details.src', sec);
    det.addEventListener('toggle', async () => {
      if (!det.open || det.dataset.loaded) return;
      det.dataset.loaded = '1';
      try {
        const r = await fetch(`../sripts/${lab.file}`);
        if (!r.ok) throw new Error(r.status);
        $('code', det).textContent = await r.text();
      } catch { $('code', det).innerHTML = `Couldn't load it here. <a href="${SRC}/${lab.file}">Open it on GitHub</a>.`; }
    });
  }
}

async function run(lab) {
  const st = state.get(lab.id);
  if (st.running) return;
  const sec = $(`#panel-${lab.id}`), btn = $('.run', sec), out = $('.out', sec);
  const inp = readInputs(lab, $('form', sec));
  if (inp.error) { out.innerHTML = `<div class="err">${esc(inp.error)}</div>`; return; }
  st.running = true; st.ran = true; st.inputs = inp;
  btn.disabled = true; btn.textContent = 'Running…';
  out.style.opacity = '.45';
  try {
    const res = await lab.run(inp);
    st.output = res.lines.join('\n');
    let body = '';
    try { body = lab.render(lab.parse(res.lines), inp); } catch (e) {
      console.warn(e);
      body = `<div class="err">The output couldn't be charted${res.truncated ? ' because it was cut off' : ''}: ${esc(e.message)}. The console output is below.</div>`;
    }
    const cut = res.truncated ? `<p class="note">The program was still printing after ${res.lines.length} lines, so it was stopped there: the iteration doesn't converge for this input.</p>` : '';
    out.innerHTML = `<div class="out-actions">${actions('results')}</div>${cut}${body}
      <div data-copyable data-title="Console output"><div class="section-h"><h2>Console output <small>· ${res.lines.length} lines, main() ran in ${res.ms < 1 ? '<1' : Math.round(res.ms)} ms</small></h2>${actions('output')}</div><pre class="console">${esc(st.output)}</pre></div>`;
    drawCharts(out);
  } catch (e) {
    console.error(e);
    st.output = '';
    out.innerHTML = `<div class="err">${esc(e.message || e)}</div>`;
  } finally {
    st.running = false; btn.disabled = false; btn.textContent = 'Run'; out.style.opacity = '';
  }
}

function select(id, push) {
  const lab = LABS.find((l) => l.id === id) || LABS[0];
  nav.querySelectorAll('button').forEach((b) => b.setAttribute('aria-selected', b.dataset.lab === lab.id));
  main.querySelectorAll('section.panel[role="tabpanel"]').forEach((s) => (s.hidden = s.id !== `panel-${lab.id}`));
  if (push) try { history.replaceState(null, '', `#${lab.id}`); } catch {}
  document.title = `Lab ${lab.tab} · Numerical Techniques Labs`;
  if (!state.get(lab.id).ran) run(lab);
  else drawCharts($(`#panel-${lab.id}`));
}

build();
select(location.hash.slice(1), false);
window.addEventListener('hashchange', () => select(location.hash.slice(1), false));
