// Orientation-independent starting positions for hydrogens: each depends only on the
// neighbouring heavy atoms, so rotating the protein cannot change which rotamer they relax into.
const XH = 0.101; // nm

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = v => { const r = Math.max(Math.hypot(...v), 1e-9); return v.map(c => c / r); };

export function placeHydrogens({ n, near, movable }, x) {
  const out = x.slice();
  const at = i => [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
  for (let parent = 0; parent < n; parent++) {
    if (movable[parent]) continue;
    const hs = near[parent].filter(i => movable[i]);
    const heavy = near[parent].filter(i => !movable[i]);
    if (!hs.length || !heavy.length) continue;
    const p = at(parent);
    const to = heavy.map(j => unit(sub(at(j), p)));
    let dirs;
    if (heavy.length >= 2) {
      const away = unit([0, 1, 2].map(k => -to.reduce((s, v) => s + v[k], 0)));
      const perp = unit(cross(to[0], to[1]));
      dirs = hs.length === 1 ? [away] : [1, -1].map(s => away.map((c, k) => 0.577 * c + s * 0.816 * perp[k]));
    } else {
      const axis = to[0];
      const beyond = near[heavy[0]].find(j => j !== parent && !movable[j]);
      const hint = beyond === undefined ? [0, 0, 0] : sub(at(beyond), at(heavy[0]));
      let side = hint.map((c, k) => c - dot(hint, axis) * axis[k]);
      if (Math.hypot(...side) < 1e-4) {
        const least = axis.map(Math.abs).indexOf(Math.min(...axis.map(Math.abs)));
        side = cross(axis, [0, 1, 2].map(k => +(k === least)));
      }
      const u = unit(side), v = cross(axis, u);
      const theta = ((hs.length === 2 ? 120 : 109.5) * Math.PI) / 180;
      dirs = ({ 1: [180], 2: [0, 180], 3: [180, 60, -60] })[hs.length].map(t => axis.map((c, k) => Math.cos(theta) * c
        + Math.sin(theta) * (Math.cos(t * Math.PI / 180) * u[k] + Math.sin(t * Math.PI / 180) * v[k])));
    }
    hs.forEach((h, i) => out.set(p.map((c, k) => c + XH * dirs[i][k]), 3 * h));
  }
  return out;
}
