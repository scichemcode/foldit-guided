// Reads the recipe's printout and rebuilds heavy-atom coordinates from the distances.
import { InputError, MAX_RESIDUES } from './topology.js';

// Rosetta heavy-atom order, which is what Foldit's band atom indices count.
export const RESIDUES = Object.fromEntries([
  ['a', 'ALA', 'CB'], ['r', 'ARG', 'CB CG CD NE CZ NH1 NH2'], ['n', 'ASN', 'CB CG OD1 ND2'],
  ['d', 'ASP', 'CB CG OD1 OD2'], ['c', 'CYS', 'CB SG'], ['q', 'GLN', 'CB CG CD OE1 NE2'],
  ['e', 'GLU', 'CB CG CD OE1 OE2'], ['g', 'GLY', ''], ['h', 'HIS', 'CB CG ND1 CD2 CE1 NE2'],
  ['i', 'ILE', 'CB CG1 CG2 CD1'], ['l', 'LEU', 'CB CG CD1 CD2'], ['k', 'LYS', 'CB CG CD CE NZ'],
  ['m', 'MET', 'CB CG SD CE'], ['f', 'PHE', 'CB CG CD1 CD2 CE1 CE2 CZ'], ['p', 'PRO', 'CB CG CD'],
  ['s', 'SER', 'CB OG'], ['t', 'THR', 'CB OG1 CG2'], ['w', 'TRP', 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2'],
  ['y', 'TYR', 'CB CG CD1 CD2 CE1 CE2 CZ OH'], ['v', 'VAL', 'CB CG1 CG2'],
].map(([code, name, side]) => [code, { name, atoms: `N CA C O ${side}`.trim().split(' ') }]));
const MAX_ANCHORS = 16;

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dist = (a, b) => Math.hypot(...sub(a, b));

function unescapeField(text) {
  try { return decodeURIComponent(text); } catch { return text; } // not UTF-8
}

function hash(text) { // cyrb53: identity only, not security
  let a = 0xdeadbeef, b = 0x41c6ce57;
  for (let i = 0; i < text.length; i++) {
    const c = text.charCodeAt(i);
    a = Math.imul(a ^ c, 2654435761);
    b = Math.imul(b ^ c, 1597334677);
  }
  a = Math.imul(a ^ (a >>> 16), 2246822507) ^ Math.imul(b ^ (b >>> 13), 3266489909);
  b = Math.imul(b ^ (b >>> 16), 2246822507) ^ Math.imul(a ^ (a >>> 13), 3266489909);
  return `${(b >>> 0).toString(36)}${(a >>> 0).toString(36)}-${text.length}`;
}

export function findCaptures(text) {
  const found = [];
  let block = null;
  for (const [line] of text.matchAll(/FEX1\|[A-Za-z0-9|.,%_-]+/g)) {
    const kind = line.split('|')[1];
    if (kind === 'BEGIN') block = [line];
    else if (block) {
      block.push(line);
      if (kind === 'END') {
        const capture = header(block);
        if (capture) found.push(capture);
        block = null;
      }
    }
  }
  return found;
}

function header(lines) {
  const f = lines[0].split('|').slice(3).map(unescapeField);
  if (f.length !== 7) return null;
  const [puzzle, name, track, move, scoreText, skipped, sym] = f;
  const score = Number(scoreText);
  if (!scoreText || !Number.isFinite(score)) return null;
  const parts = [];
  let description = '';
  for (const line of lines) {
    const r = line.split('|');
    if (r[1] === 'DESC') description = unescapeField(r[2] ?? '');
    if (r[1] === 'PART') {
      const value = Number(unescapeField(r[3] ?? ''));
      if (Number.isFinite(value)) parts.push([unescapeField(r[2]), value]);
    }
  }
  const warnings = [];
  if (Number(skipped)) warnings.push(`${Number(skipped)} non-protein segment(s) are left out of the energy.`);
  if (Number(sym)) warnings.push('Symmetry puzzle: only the base unit is analysed, without partner contacts.');
  const raw = lines.join('\n');
  return { key: hash(raw), puzzle: puzzle || name, name, track, move, score, description, parts, warnings, raw };
}

// Jacobi eigen-decomposition of a symmetric matrix, largest eigenvalue first.
function eigen(matrix) {
  const n = matrix.length, A = matrix.map(r => r.slice());
  const V = A.map((r, i) => r.map((_, j) => +(i === j)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0, all = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) { all += A[i][j] ** 2; if (i !== j) off += A[i][j] ** 2; }
    }
    if (off <= 1e-28 * all) break;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        if (!A[p][q]) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = (th < 0 ? -1 : 1) / (Math.abs(th) + Math.hypot(th, 1)), c = 1 / Math.hypot(t, 1), s = t * c;
        for (let k = 0; k < n; k++) {
          const u = A[k][p], v = A[k][q];
          A[k][p] = c * u - s * v; A[k][q] = s * u + c * v;
        }
        for (let k = 0; k < n; k++) {
          const u = A[p][k], v = A[q][k];
          A[p][k] = c * u - s * v; A[q][k] = s * u + c * v;
        }
        for (let k = 0; k < n; k++) {
          const u = V[k][p], v = V[k][q];
          V[k][p] = c * u - s * v; V[k][q] = s * u + c * v;
        }
      }
    }
  }
  return A.map((r, i) => ({ value: r[i], vector: V.map(row => row[i]) })).sort((a, b) => b.value - a.value);
}

function readDistances(text, k) {
  const d = (text ?? '').split(',').map(v => (v === '-' ? null : Number(v)));
  if (d.length !== k || d.some(v => v !== null && !(v >= 0 && v < 1e5))) throw new InputError('A distance row is damaged.');
  return d;
}

function anchorCoordinates(d) {
  const k = d.length;
  if (d.some(r => r.includes(null))) throw new InputError('Foldit refused a band between two anchors.');
  const sq = d.map((r, i) => r.map((v, j) => ((v + d[j][i]) / 2) ** 2));
  const mean = sq.map(r => r.reduce((a, b) => a + b, 0) / k), all = mean.reduce((a, b) => a + b, 0) / k;
  const top = eigen(sq.map((r, i) => r.map((v, j) => -0.5 * (v - mean[i] - mean[j] + all)))).slice(0, 3);
  if (top[2].value < 1e-2) throw new InputError('The anchor atoms are coplanar, so 3D coordinates are undetermined.');
  return d.map((_, i) => top.map(({ value, vector }) => vector[i] * Math.sqrt(value)));
}

function leastSquares(rows, rhs) {
  const M = [0, 1, 2].map(a => [0, 1, 2].map(b => rows.reduce((s, r) => s + r[a] * r[b], 0)));
  const y = [0, 1, 2].map(a => rows.reduce((s, r, i) => s + r[a] * rhs[i], 0));
  const x = [0, 0, 0];
  for (const { value, vector } of eigen(M)) {
    if (value < 1e-9) throw new InputError('An atom sees only coplanar anchors.');
    const c = (vector[0] * y[0] + vector[1] * y[1] + vector[2] * y[2]) / value;
    for (let k = 0; k < 3; k++) x[k] += c * vector[k];
  }
  return x;
}

function locate(anchors, d) {
  const use = anchors.map((_, i) => i).filter(i => d[i] !== null);
  if (use.length < 4) throw new InputError('Fewer than four distances for an atom.');
  const a = use.map(i => anchors[i]), r = use.map(i => d[i]);
  const centre = [0, 1, 2].map(k => a.reduce((s, p) => s + p[k], 0) / a.length);
  const value = a.map((p, i) => p[0] ** 2 + p[1] ** 2 + p[2] ** 2 - r[i] ** 2);
  const mean = value.reduce((s, v) => s + v, 0) / value.length;
  let x = leastSquares(a.map(p => sub(p, centre).map(v => 2 * v)), value.map(v => v - mean));
  for (let step = 0; step < 6; step++) { // Gauss-Newton on the distances themselves
    const len = a.map(p => Math.max(dist(x, p), 1e-9));
    const dx = leastSquares(a.map((p, i) => sub(x, p).map(v => v / len[i])), r.map((v, i) => v - len[i]));
    x = x.map((v, k) => v + dx[k]);
  }
  if (a.some((p, i) => Math.abs(dist(x, p) - r[i]) > 0.02)) {
    throw new InputError('Distances are inconsistent: did the protein move during capture?');
  }
  return x;
}

export function rebuild(raw) {
  const rec = raw.split('\n').map(line => line.split('|'));
  const rows = rec.filter(r => r[1] === 'ANCHOR'), atoms = rec.filter(r => r[1] === 'ATOM');
  if (rec.at(-1)[2] !== String(atoms.length)) throw new InputError('The capture is truncated.');
  const k = rows.length;
  if (k < 6 || k > MAX_ANCHORS) throw new InputError('Unexpected number of anchors.');
  if (atoms.length > MAX_RESIDUES * 15) throw new InputError(`More than ${MAX_RESIDUES} residues.`);
  const anchors = anchorCoordinates(rows.map(r => readDistances(r[3], k)));

  const segments = new Map();
  for (const [, , seg, aa, index, text] of atoms) {
    if (!Object.hasOwn(RESIDUES, aa) || !/^\d{1,4}$/.test(seg)) throw new InputError('Unknown residue in the capture.');
    if (!segments.has(seg)) segments.set(seg, { aa, xyz: new Map() });
    segments.get(seg).xyz.set(Number(index), locate(anchors, readDistances(text, k)));
  }
  if (segments.size > MAX_RESIDUES) throw new InputError(`More than ${MAX_RESIDUES} residues.`);

  let chain = 0, hand = 0;
  const residues = [...segments].sort(([a], [b]) => a - b).map(([seg, { aa, xyz }]) => {
    const { name, atoms: base } = RESIDUES[aa], n = base.length, names = base.slice();
    for (let i = 1; i <= n; i++) if (!xyz.has(i)) throw new InputError(`Segment ${seg} is missing atoms.`);
    if (xyz.has(n + 1)) { // a chain end has OXT either right after O or after the side chain
      const slot = [5, n + 1].find(s => xyz.has(s) && dist(xyz.get(s), xyz.get(3)) < 1.4);
      if (slot) names.splice(slot - 1, 0, 'OXT');
    }
    const residue = { id: seg, name, chain, atoms: names.map((atom, i) => [atom, xyz.get(i + 1)]) };
    if (names.includes('OXT')) chain++;
    const at = new Map(residue.atoms);
    if (at.has('CB')) {
      const ca = at.get('CA'), [u, v, w] = ['N', 'C', 'CB'].map(a => sub(at.get(a), ca));
      hand += Math.sign((u[1] * v[2] - u[2] * v[1]) * w[0] + (u[2] * v[0] - u[0] * v[2]) * w[1] + (u[0] * v[1] - u[1] * v[0]) * w[2]);
    }
    return residue;
  });
  if (!hand) throw new InputError('Cannot determine handedness.');
  // Distances cannot tell a structure from its mirror image; L-amino acids can.
  if (hand < 0) for (const r of residues) r.atoms = r.atoms.map(([a, [x, y, z]]) => [a, [x, y, -z]]);

  checkChemistry(residues);
  return { residues, warnings: gapWarnings(residues) };
}

// A wrong atom-index convention must fail loudly, never yield an energy.
function checkChemistry(residues) {
  const bad = residues.filter(({ atoms }) => {
    const at = new Map(atoms);
    if ([['N', 'CA', 1.46], ['CA', 'C', 1.52], ['C', 'O', 1.23]].some(([a, b, r0]) => !(Math.abs(dist(at.get(a), at.get(b)) - r0) < 0.3))) return true;
    return atoms.some(([, p], i) => i >= 4 && !atoms.slice(0, i).some(([, q]) => dist(p, q) < 2.2));
  });
  if (bad.length > Math.max(1, Math.floor(residues.length / 20))) {
    throw new InputError(`Impossible bond lengths in segments ${bad.slice(0, 8).map(r => r.id).join(', ')}. `
      + "Foldit's atom numbering may differ from what this tool assumes; please report it.");
  }
}

export function gapWarnings(residues) {
  return residues.slice(1).flatMap((next, i) => {
    const here = residues[i];
    if (here.chain !== next.chain) return [];
    const c = here.atoms.find(([a]) => a === 'C')?.[1], n = next.atoms.find(([a]) => a === 'N')?.[1];
    return c && n && dist(c, n) > 2
      ? [`Cut or gap between ${here.id} and ${next.id}: bond and angle energy there is not physical.`] : [];
  });
}
