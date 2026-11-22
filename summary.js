// One move, in plain words.
export const NOISE = 0.05; // kcal/mol; smaller changes are not judged

export function summarise(delta, scoreChange, terms) {
  const size = Math.abs(delta.total);
  const tone = size < NOISE ? 'same' : delta.total > 0 ? 'worse' : 'better';
  const headline = {
    same: 'This move hardly changed the energy',
    worse: `This move made the energy worse, by ${size.toFixed(2)} kcal/mol`,
    better: `This move made the energy better, by ${size.toFixed(2)} kcal/mol`,
  }[tone];
  const lines = [];

  const pick = compare => terms.reduce((a, b) => (compare(delta.terms[b.key], delta.terms[a.key]) ? b : a));
  const worst = pick((a, b) => a > b), best = pick((a, b) => a < b);
  if (tone === 'worse') {
    lines.push(['Main reason', `${worst.worse}${delta.where ? `, mostly at ${delta.where}` : ''}.`]);
    lines.push(['You could try', worst.tip]);
  } else if (tone === 'better') {
    lines.push(['Main gain', `${best.better}.`]);
    if (delta.terms[worst.key] >= NOISE) {
      lines.push(['What it cost', `${worst.worse}${delta.where ? `, mostly at ${delta.where}` : ''}.`]);
    }
  }

  if (scoreChange != null && Math.abs(scoreChange) >= 0.005 && tone !== 'same') {
    const up = scoreChange > 0, agree = up === (tone === 'better');
    lines.push(['Foldit’s score', `${up ? 'Went up' : 'Went down'} by ${Math.abs(scoreChange).toFixed(2)}, so it ${agree
      ? 'agrees.'
      : 'disagrees. The two measure energy differently, so neither is necessarily wrong. For ranking, Foldit’s score is the one that counts.'}`]);
  }
  return { tone, headline, lines };
}
