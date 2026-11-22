import { analyse, compare, TERMS } from './physics/analysis.js';
import { gapWarnings, rebuild } from './physics/capture.js';

self.onmessage = ({ data: { id, job, data } }) => {
  try {
    let result;
    if (job === 'terms') result = TERMS;
    else if (job === 'compare') result = compare(data.before, data.after);
    else {
      const built = data.raw ? rebuild(data.raw) : { residues: data.residues, warnings: gapWarnings(data.residues) };
      let last = 0;
      const done = analyse(built.residues, iteration => {
        if (Date.now() - last < 250) return;
        last = Date.now();
        self.postMessage({ id, progress: iteration });
      });
      result = { residues: built.residues, warnings: built.warnings, x: done.x, energy: done.energy };
    }
    self.postMessage({ id, result });
  } catch (error) {
    self.postMessage({ id, error: error.message });
  }
};
