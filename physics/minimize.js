// L-BFGS over hydrogen positions; heavy atoms never move.
import { createEvaluator } from './forces.js';

const HISTORY = 10, MAX_ITERATIONS = 500, TOLERANCE = 0.5; // kJ/mol/nm RMS; keeps the noise floor near 0.01 kcal/mol
const MAX_STEP = 0.01; // nm per coordinate

const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };

export function relaxHydrogens(sys, progress = () => {}) {
  const evaluate = createEvaluator(sys, { frozen: true });
  const free = [];
  for (let i = 0; i < sys.n; i++) if (sys.movable[i]) free.push(3 * i, 3 * i + 1, 3 * i + 2);
  const m = free.length, full = sys.x.slice(), fullGrad = new Float64Array(full.length);

  const value = (v, grad) => {
    for (let i = 0; i < m; i++) full[free[i]] = v[i];
    const e = evaluate(full, fullGrad).total;
    for (let i = 0; i < m; i++) grad[i] = fullGrad[free[i]];
    return e;
  };
  // OpenMM's criterion is an RMS over every coordinate, fixed atoms included.
  const target = TOLERANCE * Math.sqrt(3 * sys.n);

  let x = Float64Array.from(free, i => full[i]), g = new Float64Array(m), e = value(x, g);
  let memory = [], iteration = 0;
  for (; iteration < MAX_ITERATIONS; iteration++) {
    if (Math.sqrt(dot(g, g)) <= target) break;

    const p = g.map(v => -v), a = [];
    for (let k = memory.length - 1; k >= 0; k--) {
      const { s, y, rho } = memory[k];
      a[k] = rho * dot(s, p);
      for (let i = 0; i < m; i++) p[i] -= a[k] * y[i];
    }
    const last = memory[memory.length - 1];
    const scale = last ? 1 / (last.rho * dot(last.y, last.y)) : 1e-6;
    for (let i = 0; i < m; i++) p[i] *= scale;
    for (let k = 0; k < memory.length; k++) {
      const { s, y, rho } = memory[k], b = rho * dot(y, p);
      for (let i = 0; i < m; i++) p[i] += (a[k] - b) * s[i];
    }
    let downhill = dot(g, p);
    if (!(downhill < 0)) {
      memory = [];
      for (let i = 0; i < m; i++) p[i] = -1e-6 * g[i];
      downhill = dot(g, p);
    }
    let largest = 0;
    for (const v of p) largest = Math.max(largest, Math.abs(v));
    let step = Math.min(1, MAX_STEP / largest), trial, grad, next, accepted = false;
    for (let k = 0; k < 30; k++, step /= 2) {
      trial = x.map((v, i) => v + step * p[i]);
      grad = new Float64Array(m);
      next = value(trial, grad);
      if (next <= e + 1e-4 * step * downhill) { accepted = true; break; }
    }
    if (!accepted) {
      if (!memory.length) break; // no lower point along the gradient: numerical floor
      memory = [];
      continue;
    }
    const s = trial.map((v, i) => v - x[i]), y = grad.map((v, i) => v - g[i]), sy = dot(s, y);
    if (sy > 1e-12) {
      memory.push({ s, y, rho: 1 / sy });
      if (memory.length > HISTORY) memory.shift();
    }
    x = trial; g = grad; e = next;
    progress(iteration + 1);
  }
  for (let i = 0; i < m; i++) full[free[i]] = x[i];
  return { x: full, iterations: iteration };
}
