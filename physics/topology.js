// Builds the same molecular system as PDBFixer + ForceField.createSystem with ff14SB + GBn2:
// pH 7, HIE for histidine, disulfides by distance, no cutoff, no constraints.
// Algorithms follow openmm.app (see ../NOTICE.txt).
import P from './params.js';
import { placeHydrogens } from './hydrogens.js';

export const MAX_RESIDUES = 450;
export const OFFSET = 0.0195141;

export class InputError extends Error {}

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const unit = v => v.map(c => c / (Math.hypot(...v) || 1e-9));
const fits = (want, type) => want === -1 || want === type;

// residues: [{id, name, chain, atoms: [[name, [x, y, z]], ...]}], in angstrom.
// positions: the stored, already relaxed coordinates in nm (optional).
export function buildSystem(residues, positions) {
  if (!residues.length) throw new InputError('No protein residues found.');
  if (residues.length > MAX_RESIDUES) throw new InputError(`More than ${MAX_RESIDUES} residues.`);

  const input = residues.map((r, i) => {
    const first = !i || residues[i - 1].chain !== r.chain;
    const last = i === residues.length - 1 || residues[i + 1].chain !== r.chain;
    if (first && last) throw new InputError(`Residue ${r.name}${r.id} is a chain by itself.`);
    const at = new Map(r.atoms);
    if (last && !at.has('OXT') && at.has('CA') && at.has('C') && at.has('O')) {
      const away = unit(sub(sub(at.get('C'), at.get('CA')).map(c => c / 1.52),
        sub(at.get('O'), at.get('C')).map(c => c / 1.23)));
      at.set('OXT', at.get('C').map((c, k) => c + 1.25 * away[k]));
    }
    if (!last) at.delete('OXT');
    return { ...r, first, last, at };
  });

  // Topology.createDisulfideBonds: closest free partner within 3 angstrom.
  const paired = new Map();
  const cys = input.map((r, i) => (r.name === 'CYS' && r.at.has('SG') ? i : -1)).filter(i => i >= 0);
  cys.forEach((i, n) => {
    let best = 3, mate = -1;
    for (const j of cys.slice(0, n)) {
      const d = Math.hypot(...sub(input[i].at.get('SG'), input[j].at.get('SG')));
      if (d < best && !paired.has(j)) { best = d; mate = j; }
    }
    if (mate >= 0) { paired.set(i, mate); paired.set(mate, i); }
  });

  const names = [], type = [], charge = [], element = [], atomRes = [], templateIndex = [], start = [];
  const byName = [], bonds = [];
  input.forEach((r, i) => {
    if (r.first && r.last) return;
    const base = r.name === 'HIS' ? 'HIE' : paired.has(i) ? 'CYX' : r.name;
    const key = (r.first ? 'N' : r.last ? 'C' : '') + base;
    const t = Object.hasOwn(P.residues, key) ? P.residues[key] : null;
    if (!t) throw new InputError(`Unsupported residue ${r.name}${r.id}.`);
    const wanted = Object.keys(t.heavy);
    const missing = wanted.filter(a => !r.at.has(a));
    if (missing.length) throw new InputError(`${r.name}${r.id} is missing ${missing.join(', ')}.`);

    const own = new Map(), local = [];
    const add = (name, ti, p) => {
      local[ti] = names.length;
      own.set(name, names.length);
      names.push(name);
      type.push(t.atoms[ti][0]);
      charge.push(t.atoms[ti][1]);
      element.push(P.elements[t.atoms[ti][0]]);
      atomRes.push(i);
      templateIndex.push(ti);
      start.push(...p.map(c => c / 10));
    };
    // Modeller's atom order: each heavy atom is followed by its hydrogens.
    for (const name of r.at.keys()) {
      if (!Object.hasOwn(t.heavy, name)) continue;
      add(name, t.heavy[name], r.at.get(name));
      for (const [h, parent, ti] of t.hydrogens) if (parent === name) add(h, ti, r.at.get(name));
    }
    for (const [a, b] of t.bonds) bonds.push([local[a], local[b]]);
    if (!r.first) bonds.push([byName[i - 1].get('C'), own.get('N')]);
    byName.push(own);
  });
  for (const [i, j] of paired) if (i > j) bonds.push([byName[j].get('SG'), byName[i].get('SG')]);

  const n = names.length;
  const near = Array.from({ length: n }, () => []);
  for (const [a, b] of bonds) { near[a].push(b); near[b].push(a); }
  for (const list of near) list.sort((a, b) => a - b);

  const resLabels = input.map(r => `${r.name}${r.id}`);
  const sys = {
    n, names, element, atomRes, near, resLabels,
    resNames: input.map(r => r.name),
    labels: names.map((a, i) => `${resLabels[atomRes[i]]}:${a}`),
    q: Float64Array.from(charge),
    sigma: Float64Array.from(type, t => P.lj[t][0]),
    eps: Float64Array.from(type, t => P.lj[t][1]),
    movable: Uint8Array.from(element, e => (e === 'H' ? 1 : 0)),
  };
  sys.signature = JSON.stringify([sys.labels, bonds.map(([a, b]) => (a < b ? [a, b] : [b, a])).sort()]);

  sys.bonds = bonds.map(([a, b]) => {
    const d = P.bonds.find(([u, v]) => (u === type[a] && v === type[b]) || (v === type[a] && u === type[b]));
    if (!d) throw new Error(`No bond parameters for ${sys.labels[a]}-${sys.labels[b]}`);
    return [a, b, d[2], d[3]];
  });

  const angles = [];
  for (let j = 0; j < n; j++) {
    for (let a = 0; a < near[j].length; a++) {
      for (let b = a + 1; b < near[j].length; b++) angles.push([near[j][a], j, near[j][b]]);
    }
  }
  sys.angles = angles.map(([i, j, k]) => {
    const [a, b, c] = [type[i], type[j], type[k]];
    const d = P.angles.find(([u, v, w]) => v === b && ((u === a && w === c) || (u === c && w === a)));
    if (!d) throw new Error(`No angle parameters for ${sys.labels[i]}-${sys.labels[j]}-${sys.labels[k]}`);
    return [i, j, k, d[3], d[4]];
  });

  sys.torsions = [...properTorsions(type, near, bonds), ...improperTorsions(type, near, element, atomRes, templateIndex)];
  addExclusions(sys);
  addGeneralizedBorn(sys);

  if (positions) {
    if (positions.length !== 3 * n) throw new InputError('Stored coordinates do not fit this protein.');
    sys.x = Float64Array.from(positions);
  } else {
    sys.x = placeHydrogens(sys, Float64Array.from(start));
  }
  return sys;
}

function properTorsions(type, near, bonds) {
  const out = [], cache = new Map();
  const find = (a, b, c, d) => {
    let match = null;
    for (const def of P.propers) {
      const [u, v, w, z] = def;
      const forward = fits(v, b) && fits(w, c) && fits(z, d) && fits(u, a);
      if (!forward && !(fits(w, b) && fits(v, c) && fits(u, d) && fits(z, a))) continue;
      const wild = u < 0 || v < 0 || w < 0 || z < 0;
      if (!match || !wild) match = def; // specific beats wildcard
      if (!wild) break;
    }
    return match;
  };
  for (const [b, c] of bonds) {
    for (const a of near[b]) {
      for (const d of near[c]) {
        if (a === c || d === b || a === d) continue;
        const atoms = a < d ? [a, b, c, d] : [d, c, b, a];
        const key = atoms.map(i => type[i]).join();
        if (!cache.has(key)) cache.set(key, find(...atoms.map(i => type[i])));
        for (const [m, phase, k] of cache.get(key)?.[4] ?? []) if (k) out.push([...atoms, m, phase, k]);
      }
    }
  }
  return out;
}

// forcefield._matchImproper with ordering="amber", including its per-type-signature cache.
const PERMUTATIONS = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]];
function improperTorsions(type, near, element, res, ti) {
  const out = [], cache = new Map();
  const after = (a, b) => res[a] > res[b] || (res[a] === res[b] && ti[a] > ti[b]);

  const match = atoms => {
    let found = null;
    for (const def of P.impropers) {
      const [u, v, w, z] = def;
      const wild = u < 0 || v < 0 || w < 0 || z < 0;
      if (found && wild) continue;
      if (!fits(u, type[atoms[0]])) continue;
      for (const [p2, p3, p4] of PERMUTATIONS) {
        let a2 = atoms[p2], a3 = atoms[p3], a4 = atoms[p4];
        const t2 = type[a2], t3 = type[a3], t4 = type[a4];
        if (!fits(v, t2) || !fits(w, t3) || !fits(z, t4)) continue;
        if (!wild) {
          if (t2 === t4 && after(a2, a4)) [a2, a4] = [a4, a2];
          if (t3 === t4 && after(a3, a4)) [a3, a4] = [a4, a3];
          if (t2 === t3 && after(a2, a3)) [a2, a3] = [a3, a2];
        } else {
          const e2 = element[a2], e3 = element[a3], e4 = element[a4];
          if (e2 === e4 && after(a2, a4)) [a2, a4] = [a4, a2];
          if (e3 === e4 && after(a3, a4)) [a3, a4] = [a4, a3];
          if (after(a2, a3)) [a2, a3] = [a3, a2];
        }
        found = { order: [a2, a3, atoms[0], a4].map(a => atoms.indexOf(a)), terms: def[4] };
        break;
      }
    }
    return found;
  };

  near.forEach((list, centre) => {
    for (let a = 0; a < list.length; a++) {
      for (let b = a + 1; b < list.length; b++) {
        for (let c = b + 1; c < list.length; c++) {
          const atoms = [centre, list[a], list[b], list[c]];
          const key = atoms.map(i => type[i]).join();
          if (!cache.has(key)) cache.set(key, match(atoms));
          const hit = cache.get(key);
          if (!hit) continue;
          for (const [m, phase, k] of hit.terms) if (k) out.push([...hit.order.map(p => atoms[p]), m, phase, k]);
        }
      }
    }
  });
  return out;
}

function addExclusions(sys) {
  const { n, near, q, sigma, eps } = sys;
  const [scaleQ, scaleLJ] = P.scale14;
  sys.excluded = [];
  sys.pairs14 = [];
  for (let i = 0; i < n; i++) {
    const depth = new Map([[i, 0]]);
    let edge = [i];
    for (let d = 1; d <= 3; d++) {
      const next = [];
      for (const a of edge) for (const b of near[a]) if (!depth.has(b)) { depth.set(b, d); next.push(b); }
      edge = next;
    }
    const later = [];
    for (const [j, d] of depth) {
      if (j <= i) continue;
      later.push(j);
      if (d === 3) sys.pairs14.push([i, j, q[i] * q[j] * scaleQ, (sigma[i] + sigma[j]) / 2, Math.sqrt(eps[i] * eps[j]) * scaleLJ]);
    }
    sys.excluded.push(later);
  }
}

function addGeneralizedBorn(sys) {
  const { n, element, near, names, atomRes, resNames } = sys;
  const radius = i => { // mbondi3, angstrom
    const e = element[i], name = names[i];
    if (resNames[atomRes[i]] === 'ARG' && (name.startsWith('HH') || name.startsWith('HE'))) return 1.17;
    if (e === 'H') return element[near[i][0]] === 'N' ? 1.3 : 1.2;
    if (e === 'O') {
      const c = near[i][0];
      const carboxylate = near[i].length === 1 && element[c] === 'C' && near[c].length === 3
        && near[c].some(a => a !== i && element[a] === 'O');
      return carboxylate ? 1.4 : 1.5;
    }
    return { C: 1.7, N: 1.55, S: 1.8 }[e];
  };
  const r = Array.from({ length: n }, (_, i) => radius(i) / 10);
  const unique = [...new Set(r)].sort((a, b) => a - b);
  const size = unique.length;
  const position = unique.map(v => (v - 0.1) * 200);
  const interpolate = (table, row, col) => {
    const split = p => {
      if (p <= 0) return [[0, 1]];
      if (p >= 20) return [[20, 1]];
      const f = Math.floor(p);
      return [[f, f + 1 - p], [f + 1, p - f]];
    };
    let v = 0;
    for (const [a, wa] of split(row)) for (const [b, wb] of split(col)) if (wa && wb) v += wa * wb * table[a * 21 + b];
    return v;
  };
  // Same layout as OpenMM's table; forces.js reads it column-major, like Discrete2DFunction.
  const table = full => Float64Array.from({ length: size * size },
    (_, k) => interpolate(full, position[Math.floor(k / size)], position[k % size]));
  const per = k => Float64Array.from(element, e => P.gb.elements[e][k]);
  const born = Float64Array.from(r, v => v - OFFSET);
  sys.gb = {
    size, born,
    radius: Float64Array.from(r),
    scaled: per(0).map((s, i) => s * born[i]),
    alpha: per(1), beta: per(2), gamma: per(3),
    index: Int32Array.from(r, v => unique.indexOf(v)),
    m0: table(P.gb.m0),
    d0: table(P.gb.d0),
  };
}
