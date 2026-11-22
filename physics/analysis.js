// The analysis of foldit_energy.py: energy per term, and why the dominant penalty went up.
import { buildSystem } from './topology.js';
import { COULOMB, angleOf, createEvaluator, dihedral } from './forces.js';
import { relaxHydrogens } from './minimize.js';
import { residueSasa } from './sasa.js';

export const TERMS = [
  { key: 'bond', plain: 'Stretched bonds', worse: 'Bonds were stretched or squeezed', better: 'Bond lengths are closer to normal',
    tip: 'Run Idealize on that area, or loosen any bands pulling on it.',
    label: 'Bond stretch', equation: '½ k_b (r − r₀)²', foldit: 'Ideality', hint: "'Ideality': stretched bond, usually from pulling with bands." },
  { key: 'angle', plain: 'Bent bond angles', worse: 'Bond angles were bent out of shape', better: 'Bond angles are closer to normal',
    tip: 'Wiggle that area, or loosen any bands pulling on it.',
    label: 'Angle bend', equation: '½ k_θ (θ − θ₀)²', foldit: 'Ideality', hint: "'Ideality': distorted geometry, often after a forced Rebuild or band pull." },
  { key: 'torsion', plain: 'Awkward twists', worse: 'The chain or a side chain was twisted awkwardly', better: 'Twists in the chain are more comfortable',
    tip: 'Wiggle the backbone there, or Shake the side chains.',
    label: 'Torsion', equation: 'Σ k_φ [1 + cos(nφ − δ)]', foldit: 'Backbone / Sidechain', hint: "'Backbone' / 'Sidechain': strained dihedral. Wiggle Backbone or Shake." },
  { key: 'lj_rep', plain: 'Clashes', worse: 'Atoms were pushed too close together', better: 'Clashes between atoms were relieved',
    tip: 'Shake, then Wiggle, around that spot.',
    label: 'LJ repulsion', equation: '4ε (σ/r)¹²', foldit: 'Clashing', hint: "'Clashing': try Shake or Wiggle Sidechains." },
  { key: 'lj_disp', plain: 'Packing', worse: 'Packing became looser', better: 'The protein packs more tightly',
    tip: 'Bring those parts closer, then Shake and Wiggle to fill the gap.',
    label: 'LJ dispersion', equation: '−4ε (σ/r)⁶', foldit: 'Packing', hint: "'Packing': the core may have voids. Wiggle or Shake to re-pack." },
  { key: 'coulomb', plain: 'Attraction of charges', worse: 'Favourable attractions between charges, such as hydrogen bonds, got weaker', better: 'Attractions between charges, such as hydrogen bonds, got stronger',
    tip: 'Look for a broken hydrogen bond or salt bridge there and restore it.',
    label: 'Coulomb', equation: 'q_i q_j / (4πε₀ r)', foldit: 'Bonding', hint: "'Bonding': look for lost H-bonds or broken salt bridges." },
  { key: 'gb_solv', plain: 'Contact with water', worse: 'Contact with the surrounding water became less favourable', better: 'Contact with the surrounding water became more favourable',
    tip: 'Turn greasy (orange) side chains inward and charged (blue) ones outward.',
    label: 'GB solvation', equation: 'ΔG_pol (GBn2) + ΔG_np (SA)', foldit: 'Hiding', hint: "'Hiding': hydrophobics exposed or charges buried." },
];
const KCAL = 1 / 4.184;
const CHARGE = { ASP: -1, GLU: -1, LYS: 1, ARG: 1 };
const NONPOLAR = new Set(['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'CYS', 'TYR']);
const kc = (v, d = 2) => (v * KCAL).toFixed(d);
const signed = (v, d = 2) => (v > 0 ? '+' : '') + v.toFixed(d);

export function analyse(residues, progress) {
  const sys = buildSystem(residues);
  const { x, iterations } = relaxHydrogens(sys, progress);
  const E = createEvaluator(sys)(x);
  if (!Number.isFinite(E.total)) throw new Error('The energy is not finite (overlapping atoms?).');
  const terms = Object.fromEntries(TERMS.map(t => [t.key, E[t.key] * KCAL]));
  return { x, iterations, energy: { total: E.total * KCAL, terms } };
}

export function compare(before, after) {
  const A = buildSystem(before.residues, before.x), B = buildSystem(after.residues, after.x);
  if (A.signature !== B.signature) return { comparable: false };

  const terms = Object.fromEntries(TERMS.map(({ key }) => [key, after.energy.terms[key] - before.energy.terms[key]]));
  const worst = TERMS.reduce((a, b) => (terms[b.key] > terms[a.key] ? b : a));
  const best = TERMS.reduce((a, b) => (terms[b.key] < terms[a.key] ? b : a));
  const total = after.energy.total - before.energy.total;
  const out = { comparable: true, total, terms, causes: [] };
  if (terms[best.key] < -0.005) out.gain = `Largest gain: ${best.label}, ${signed(terms[best.key])} kcal/mol.`;
  if (terms[worst.key] < 0.005) return out;

  out.heading = total > 0 ? `Why ${worst.label} went up` : `Price paid for the net gain: ${worst.label} went up`;
  out.equation = worst.equation;
  out.hint = `In Foldit this is closest to ${worst.hint}`;
  if (worst.key === 'gb_solv') out.causes = solvationCauses(A, B);
  else if (['bond', 'angle', 'torsion'].includes(worst.key)) out.causes = bondedCauses(A, B, worst.key);
  else out.causes = pairCauses(A, B, worst.key);
  if (worst.key !== 'gb_solv' && out.causes.length) out.where = out.causes[0].where;
  return out;
}

// kJ/mol per atom pair, given as [i, j, r, qq, sigma, eps]; exclusions and 1-4 scaling as in OpenMM.
function forEachPair(sys, visit, wanted = () => true) {
  const { n, x, q, sigma, eps, atomRes, excluded, pairs14 } = sys;
  const distance = (i, j) => Math.max(Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]), 1e-4);
  const skip = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    for (const j of excluded[i]) skip[j] = i;
    for (let j = i + 1; j < n; j++) {
      if (skip[j] === i || !wanted(atomRes[i], atomRes[j])) continue;
      visit(i, j, distance(i, j), q[i] * q[j], (sigma[i] + sigma[j]) / 2, Math.sqrt(eps[i] * eps[j]));
    }
  }
  for (const [i, j, qq, s, e] of pairs14) if (wanted(atomRes[i], atomRes[j])) visit(i, j, distance(i, j), qq, s, e);
}

const PAIR_ENERGY = {
  lj_rep: (r, qq, s, e) => 4 * e * (s / r) ** 12,
  lj_disp: (r, qq, s, e) => -4 * e * (s / r) ** 6,
  coulomb: (r, qq) => (COULOMB * qq) / r,
};

function pairCauses(A, B, key) {
  const R = B.resLabels.length, change = new Float64Array(R * R), energy = PAIR_ENERGY[key];
  for (const [sys, sign] of [[A, -1], [B, 1]]) {
    forEachPair(sys, (i, j, r, qq, s, e) => {
      const a = sys.atomRes[i], b = sys.atomRes[j];
      change[Math.min(a, b) * R + Math.max(a, b)] += sign * energy(r, qq, s, e);
    });
  }
  const top = [];
  change.forEach((v, k) => {
    if (v <= 0 || (top.length === 3 && v <= top[2][0])) return;
    top.push([v, k]);
    top.sort((a, b) => b[0] - a[0]);
    top.length = Math.min(top.length, 3);
  });
  return top.map(([v, k]) => {
    const ri = Math.floor(k / R), rj = k % R;
    const title = ri === rj ? B.resLabels[ri] : `${B.resLabels[ri]}–${B.resLabels[rj]}`;
    return { where: title, title: `${title}: ${signed(v * KCAL)} kcal/mol`, notes: explainPair(A, B, key, ri, rj) };
  });
}

function explainPair(A, B, key, ri, rj) {
  const wanted = (a, b) => (a === ri && b === rj) || (a === rj && b === ri);
  const old = new Map();
  forEachPair(A, (i, j, r, qq, s, e) => old.set(`${i},${j}`, { r, rep: PAIR_ENERGY.lj_rep(r, qq, s, e) }), wanted);
  let best = null, nearA = Infinity, nearB = Infinity;
  for (const { r } of old.values()) nearA = Math.min(nearA, r);
  forEachPair(B, (i, j, r, qq, s, e) => {
    const rep = PAIR_ENERGY.lj_rep(r, qq, s, e), was = old.get(`${i},${j}`);
    const rise = rep - (was?.rep ?? 0);
    nearB = Math.min(nearB, r);
    if (!best || rise > best.rise) best = { i, j, r, s, e, rep, was, rise };
  }, wanted);

  if (key === 'lj_rep') {
    const { i, j, r, s, e, rep, was } = best, ra = r * 10;
    const notes = [
      `${B.labels[i]} ··· ${B.labels[j]}: r = ${ra.toFixed(2)} Å (σ = ${(s * 10).toFixed(2)} Å, r_min = 2^(1/6)σ = ${(2 ** (1 / 6) * s * 10).toFixed(2)} Å)`,
      `r/σ = ${(r / s).toFixed(2)}, (σ/r)¹² = ${((s / r) ** 12).toPrecision(3)}, E_rep = 4ε(σ/r)¹² = ${kc(rep)} kcal/mol (ε = ${kc(e, 3)} kcal/mol)`,
    ];
    if (was) notes.push(`Before the move r = ${(was.r * 10).toFixed(2)} Å (E_rep ${kc(was.rep)} → ${kc(rep)} kcal/mol).`);
    if (ra > 0.2) notes.push(`Pauli exchange repulsion. On the r⁻¹² wall, 0.1 Å closer would multiply E_rep by ${((ra / (ra - 0.1)) ** 12).toFixed(2)}×.`);
    return notes;
  }

  const contact = `Closest contact ${(nearA * 10).toFixed(2)} Å → ${(nearB * 10).toFixed(2)} Å.`;
  if (key === 'lj_disp') {
    return [contact, nearB > nearA
      ? 'London dispersion (−C₆/r⁶) lost: the two groups moved apart, leaving looser packing.'
      : 'London dispersion (−C₆/r⁶) lost: fewer atom pairs remain in good contact.'];
  }
  const sign = (CHARGE[B.resNames[ri]] ?? 0) * (CHARGE[B.resNames[rj]] ?? 0);
  const what = ri !== rj && sign < 0 ? 'A salt bridge weakened (opposite charges).'
    : ri !== rj && sign > 0 ? 'Like charges became less favourable.'
    : 'Partial charges rearranged (for example an H-bond broke).';
  return [contact, what, 'Coulomb here is unscreened (ε = 1); GB partly screens it, so read the two terms together.'];
}

function bondedCauses(A, B, key) {
  const atoms = t => t.map(i => B.names[i]).join('-');
  let rows;
  if (key === 'bond') {
    const length = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]);
    rows = B.bonds.map(([i, j, r0, k]) => {
      const a = length(A.x, i, j), b = length(B.x, i, j);
      return { t: [i, j], rise: 0.5 * k * ((b - r0) ** 2 - (a - r0) ** 2), note: `r ${(a * 10).toFixed(3)} → ${(b * 10).toFixed(3)} Å (r₀ ${(r0 * 10).toFixed(3)} Å)` };
    });
  } else if (key === 'angle') {
    const deg = v => ((v * 180) / Math.PI).toFixed(1);
    rows = B.angles.map(([i, j, l, t0, k]) => {
      const a = angleOf(A.x, i, j, l), b = angleOf(B.x, i, j, l);
      return { t: [i, j, l], rise: 0.5 * k * ((b - t0) ** 2 - (a - t0) ** 2), note: `θ ${deg(a)}° → ${deg(b)}° (θ₀ ${deg(t0)}°)` };
    });
  } else {
    const grouped = new Map(), deg = v => ((v * 180) / Math.PI).toFixed(0);
    for (const [i, j, l, m, n, phase, k] of B.torsions) {
      const id = `${i},${j},${l},${m}`;
      if (!grouped.has(id)) grouped.set(id, { t: [i, j, l, m], a: dihedral(A.x, i, j, l, m), b: dihedral(B.x, i, j, l, m), series: [], rise: 0 });
      const g = grouped.get(id);
      g.series.push(`n=${n}, δ=${deg(phase)}°, k=${kc(k)}`);
      g.rise += k * (Math.cos(n * g.b - phase) - Math.cos(n * g.a - phase));
    }
    rows = [...grouped.values()].map(g => ({ ...g, note: `φ ${deg(g.a)}° → ${deg(g.b)}° (${g.series.join('; ')})` }));
  }
  return rows.filter(r => r.rise > 0).sort((a, b) => b.rise - a.rise).slice(0, 3).map(({ t, rise, note }) =>
    ({ where: B.resLabels[B.atomRes[t[1]]], title: `${B.resLabels[B.atomRes[t[1]]]} ${atoms(t)}: ${signed(rise * KCAL)} kcal/mol`, notes: [note] }));
}

function solvationCauses(A, B) {
  const a = residueSasa(A), b = residueSasa(B);
  const rows = B.resLabels.map((label, i) => ({ label, name: B.resNames[i], d: b[i] - a[i] }));
  const show = list => list.slice(0, 3).map(r => `${r.label} (${signed(r.d, 0)} Å²)`);
  const exposed = rows.filter(r => NONPOLAR.has(r.name) && r.d > 5).sort((p, q) => q.d - p.d);
  const buried = rows.filter(r => r.name in CHARGE && r.d < -5).sort((p, q) => p.d - q.d);
  const causes = [];
  if (exposed.length) causes.push({ title: 'Nonpolar residues newly exposed to water', notes: show(exposed) });
  if (buried.length) causes.push({ title: 'Charged residues newly buried (desolvation penalty)', notes: show(buried) });
  causes.push({ title: 'What this term is', notes: ['Generalized Born screening by continuum water (ε = 78.5) plus a surface-area term for the hydrophobic effect.'] });
  return causes;
}
