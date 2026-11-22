// Export and import of a player's history. An imported file is untrusted.
const FORMAT = 'foldit-move-energetics/1';

const text = (v, max = 200) => (typeof v === 'string' ? v.slice(0, max) : '');
const number = v => (typeof v === 'number' && Number.isFinite(v) ? v : null);

export function exportSteps(steps) {
  return JSON.stringify({ format: FORMAT, steps: steps.map(s => ({ ...s, x: s.x && Array.from(s.x) })) });
}

export function importSteps(source) {
  const file = JSON.parse(source);
  if (file?.format !== FORMAT || !Array.isArray(file.steps)) throw new Error('This is not an exported history file.');
  return file.steps.map(s => {
    const residues = Array.isArray(s.residues) ? s.residues.map(r => {
      const atoms = r.atoms.map(([name, p]) => {
        if (!Array.isArray(p) || p.length !== 3 || p.some(v => number(v) === null)) throw new Error('The file has damaged coordinates.');
        return [text(name, 4), p.slice()];
      });
      return { id: text(String(r.id), 12), name: text(r.name, 3), chain: number(r.chain) ?? 0, atoms };
    }) : undefined;
    const raw = text(s.raw, 5e6) || undefined;
    if (!text(s.key) || !text(s.challenge) || (!residues && !raw)) throw new Error('The file has a damaged step.');
    // Energies are recomputed, never trusted from the file.
    return {
      key: text(s.key), challenge: text(s.challenge), name: text(s.name), track: text(s.track), move: text(s.move),
      description: text(s.description, 5000), order: number(s.order) ?? Date.now(), score: number(s.score),
      parts: (Array.isArray(s.parts) ? s.parts : []).filter(p => number(p?.[1]) !== null).map(([k, v]) => [text(k, 60), v]),
      notes: (Array.isArray(s.notes) ? s.notes : []).map(n => text(n, 400)),
      status: 'pending', residues, raw,
    };
  });
}
