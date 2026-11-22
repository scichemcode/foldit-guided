import { RESIDUES } from './capture.js';
import { InputError } from './topology.js';

const KNOWN = new Set(Object.values(RESIDUES).map(r => r.name));
const ALIAS = { HID: 'HIS', HIE: 'HIS', HIP: 'HIS', HSD: 'HIS', HSE: 'HIS', CYX: 'CYS' };

export function parsePdb(text) {
  const list = [], skipped = new Set();
  let ter = 0, last = '', current = null, het = false;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('ENDMDL')) break;
    if (line.startsWith('TER')) { ter++; last = ''; continue; }
    if (line.startsWith('HETATM')) { het = true; continue; }
    if (!line.startsWith('ATOM  ')) continue;
    const element = line.slice(76, 78).trim().toUpperCase();
    let atom = line.slice(12, 16).trim();
    if (element ? element === 'H' || element === 'D' : /^\d*H/.test(atom)) continue;
    const code = line.slice(17, 20).trim(), name = ALIAS[code] ?? code;
    if (!KNOWN.has(name)) { skipped.add(code); continue; }
    const p = [30, 38, 46].map(i => parseFloat(line.slice(i, i + 8)));
    if (!p.every(Number.isFinite)) throw new InputError('An ATOM record has unreadable coordinates.');
    const key = ter + line.slice(21, 27);
    if (key !== last) {
      current = { name, seq: line.slice(22, 27).trim(), letter: line[21].trim(), ter, at: new Map() };
      list.push(current);
      last = key;
    }
    if (name === 'ILE' && atom === 'CD') atom = 'CD1';
    if (atom === 'OT1') atom = 'O';
    if (atom === 'OT2') atom = 'OXT';
    if (!current.at.has(atom)) current.at.set(atom, p); // first alternate location wins
  }
  if (!list.length) {
    throw new InputError("No protein ATOM records found. Foldit's .ir_solution saves are binary: use the recipe for those.");
  }

  const warnings = [];
  if (het) warnings.push('Ligands, waters and other HETATM records are left out.');
  if (skipped.size) warnings.push(`Left out unsupported residues: ${[...skipped].slice(0, 8).join(', ')}.`);
  const many = new Set(list.map(r => r.letter)).size > 1;
  let chain = 0;
  const residues = list.map((r, i) => {
    const before = list[i - 1];
    if (before) {
      const c = before.at.get('C'), n = r.at.get('N');
      const gap = c && n && Math.hypot(...c.map((v, k) => v - n[k])) > 2;
      if (gap && before.ter === r.ter && before.letter === r.letter) {
        warnings.push(`Gap between ${before.seq} and ${r.seq}: treated as two chains.`);
      }
      if (gap || before.ter !== r.ter || before.letter !== r.letter) chain++;
    }
    return { id: many ? `${r.seq}${r.letter}` : r.seq, name: r.name, chain, atoms: [...r.at] };
  });
  return { residues, warnings };
}
