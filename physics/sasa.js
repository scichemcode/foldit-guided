// Shrake-Rupley surface area per residue, as Biopython computes it.
const RADIUS = { H: 1.2, C: 1.7, N: 1.55, O: 1.52, S: 1.8 }, PROBE = 1.4, POINTS = 100;

const SPHERE = Array.from({ length: POINTS }, (_, k) => {
  const z = 1 - (2 * k + 1) / POINTS, r = Math.sqrt(1 - z * z), t = k * Math.PI * (3 - Math.sqrt(5));
  return [Math.cos(t) * r, Math.sin(t) * r, z];
});

export function residueSasa({ n, element, atomRes, resLabels, x }) {
  const total = new Float64Array(resLabels.length);
  const R = element.map(e => RADIUS[e] + PROBE);
  for (let i = 0; i < n; i++) {
    const near = [];
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const dx = 10 * (x[3 * j] - x[3 * i]), dy = 10 * (x[3 * j + 1] - x[3 * i + 1]), dz = 10 * (x[3 * j + 2] - x[3 * i + 2]);
      if (dx * dx + dy * dy + dz * dz < (R[i] + R[j]) ** 2) near.push([dx, dy, dz, R[j] * R[j]]);
    }
    let open = 0, hit = 0;
    for (const [sx, sy, sz] of SPHERE) {
      const px = sx * R[i], py = sy * R[i], pz = sz * R[i];
      let buried = false;
      for (let t = 0; t < near.length; t++) {
        const k = (hit + t) % near.length, o = near[k]; // the last occluder usually hides the next point too
        if ((px - o[0]) ** 2 + (py - o[1]) ** 2 + (pz - o[2]) ** 2 < o[3]) { buried = true; hit = k; break; }
      }
      if (!buried) open++;
    }
    total[atomRes[i]] += (open / POINTS) * 4 * Math.PI * R[i] * R[i];
  }
  return total;
}
