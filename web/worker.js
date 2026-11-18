// Runs one lab's compiled main() in a fresh WebAssembly instance and returns what it printed.
// The page passes optional values for the lab's mutable globals (float M[4][5] in 2-lab,
// float X[5] and Y[5] in 5-lab); they are written into linear memory before main() starts,
// exactly as if the arrays had been initialised with them.
const LINE_LIMIT = 6000;

self.onmessage = async ({ data: { lab, globals } }) => {
  const lines = [];
  let truncated = false;
  const print = (s) => {
    lines.push(s);
    if (lines.length >= LINE_LIMIT) { truncated = true; throw new Error('output limit'); }
  };
  try {
    const { default: createLab } = await import(`./wasm/lab${lab}.mjs`);
    const mod = await createLab({ print, printErr: (s) => print(`[stderr] ${s}`) });
    for (const [name, values] of Object.entries(globals || {})) {
      const addr = mod[`_${name}`];
      if (typeof addr !== 'number') throw new Error(`lab ${lab} has no exported global ${name}`);
      mod.HEAPF32.set(Float32Array.from(values), addr >> 2);
    }
    const t0 = performance.now();
    let code = 0;
    try {
      code = mod.callMain([]);
    } catch (e) {
      if (!truncated) throw e;
    }
    const ms = performance.now() - t0;
    self.postMessage({ ok: true, lines, truncated, code, ms });
  } catch (e) {
    self.postMessage({ ok: false, error: String(e && e.message || e), lines });
  }
};
