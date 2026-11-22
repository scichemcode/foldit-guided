// Energy (kJ/mol) and gradient (kJ/mol/nm), matching OpenMM's NoCutoff forces.
import { OFFSET } from './topology.js';

export const COULOMB = 138.93545764438198;
const GB_K = 138.935485 * (1 - 1 / 78.5);
const NECK_SCALE = 0.826836, NECK_CUT = 0.68;

export function angleOf(x, i, j, k) {
  const a = 3 * i, b = 3 * j, c = 3 * k;
  const ux = x[a] - x[b], uy = x[a + 1] - x[b + 1], uz = x[a + 2] - x[b + 2];
  const vx = x[c] - x[b], vy = x[c + 1] - x[b + 1], vz = x[c + 2] - x[b + 2];
  return Math.atan2(Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx), ux * vx + uy * vy + uz * vz);
}

// IUPAC sign. With a gradient array g, also adds dE/dphi(phi) times dphi/dx.
export function dihedral(x, i, j, k, l, g, slope) {
  const a = 3 * i, b = 3 * j, c = 3 * k, d = 3 * l;
  const fx = x[a] - x[b], fy = x[a + 1] - x[b + 1], fz = x[a + 2] - x[b + 2]; // ij
  const gx = x[c] - x[b], gy = x[c + 1] - x[b + 1], gz = x[c + 2] - x[b + 2]; // kj
  const hx = x[c] - x[d], hy = x[c + 1] - x[d + 1], hz = x[c + 2] - x[d + 2]; // kl
  const mx = fy * gz - fz * gy, my = fz * gx - fx * gz, mz = fx * gy - fy * gx;
  const nx = gy * hz - gz * hy, ny = gz * hx - gx * hz, nz = gx * hy - gy * hx;
  const len2 = gx * gx + gy * gy + gz * gz, len = Math.sqrt(len2);
  const phi = Math.atan2(len * (fx * nx + fy * ny + fz * nz), mx * nx + my * ny + mz * nz);
  if (g) {
    const m2 = mx * mx + my * my + mz * mz, n2 = nx * nx + ny * ny + nz * nz;
    if (m2 < 1e-20 || n2 < 1e-20) return phi; // collinear: direction undefined
    const s = slope(phi), si = (s * len) / m2, sl = (-s * len) / n2;
    const p = (fx * gx + fy * gy + fz * gz) / len2, q = (hx * gx + hy * gy + hz * gz) / len2;
    const ix = si * mx, iy = si * my, iz = si * mz, lx = sl * nx, ly = sl * ny, lz = sl * nz;
    g[a] += ix; g[a + 1] += iy; g[a + 2] += iz;
    g[d] += lx; g[d + 1] += ly; g[d + 2] += lz;
    g[b] += (p - 1) * ix - q * lx; g[b + 1] += (p - 1) * iy - q * ly; g[b + 2] += (p - 1) * iz - q * lz;
    g[c] += (q - 1) * lx - p * ix; g[c + 1] += (q - 1) * ly - p * iy; g[c + 2] += (q - 1) * lz - p * iz;
  }
  return phi;
}

function bondedForces(sys, x, g, E) {
  for (const [i, j, r0, k] of sys.bonds) {
    const a = 3 * i, b = 3 * j;
    const dx = x[a] - x[b], dy = x[a + 1] - x[b + 1], dz = x[a + 2] - x[b + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz), f = (k * (r - r0)) / r;
    E.bond += 0.5 * k * (r - r0) * (r - r0);
    g[a] += f * dx; g[a + 1] += f * dy; g[a + 2] += f * dz;
    g[b] -= f * dx; g[b + 1] -= f * dy; g[b + 2] -= f * dz;
  }

  for (const [i, j, k, t0, kt] of sys.angles) {
    const a = 3 * i, b = 3 * j, c = 3 * k;
    const ux = x[a] - x[b], uy = x[a + 1] - x[b + 1], uz = x[a + 2] - x[b + 2];
    const vx = x[c] - x[b], vy = x[c + 1] - x[b + 1], vz = x[c + 2] - x[b + 2];
    const px = uy * vz - uz * vy, py = uz * vx - ux * vz, pz = ux * vy - uy * vx;
    const p = Math.max(Math.sqrt(px * px + py * py + pz * pz), 1e-8);
    const t = Math.atan2(p, ux * vx + uy * vy + uz * vz), d = kt * (t - t0);
    E.angle += 0.5 * d * (t - t0);
    const su = d / ((ux * ux + uy * uy + uz * uz) * p), sv = -d / ((vx * vx + vy * vy + vz * vz) * p);
    const ax = su * (uy * pz - uz * py), ay = su * (uz * px - ux * pz), az = su * (ux * py - uy * px);
    const cx = sv * (vy * pz - vz * py), cy = sv * (vz * px - vx * pz), cz = sv * (vx * py - vy * px);
    g[a] += ax; g[a + 1] += ay; g[a + 2] += az;
    g[c] += cx; g[c + 1] += cy; g[c + 2] += cz;
    g[b] -= ax + cx; g[b + 1] -= ay + cy; g[b + 2] -= az + cz;
  }

  let m, phase, k;
  const slope = phi => -k * m * Math.sin(m * phi - phase);
  for (const t of sys.torsions) {
    m = t[4]; phase = t[5]; k = t[6];
    const phi = dihedral(x, t[0], t[1], t[2], t[3], g, slope);
    E.torsion += k * (1 + Math.cos(m * phi - phase));
  }
}

// With frozen = true only gradients on hydrogens are valid; pairs of fixed atoms are then cheaper.
export function createEvaluator(sys, { frozen = false } = {}) {
  const { n, q, sigma, eps, excluded, pairs14 } = sys;
  const { born, radius, scaled, alpha, beta, gamma, index, size, m0, d0 } = sys.gb;
  const moving = frozen ? sys.movable : new Uint8Array(n).fill(1);
  const halfSigma = sigma.map(s => s / 2), rootEps = eps.map(Math.sqrt);

  let count = 0;
  for (let i = 0; i < n; i++) if (!moving[i]) count++;
  const livePairs = (n * (n - 1) - count * (count - 1)) / 2;
  const Store = n > 2500 ? Float32Array : Float64Array; // memory, for big proteins
  const dI = new Store(2 * livePairs);

  const skip = new Int32Array(n).fill(-1);
  const I = new Float64Array(n), B = new Float64Array(n), dBdI = new Float64Array(n), dEdB = new Float64Array(n);
  let fixedI = null, fixedRep = 0, fixedDisp = 0, fixedCoul = 0;
  let slope = 0; // second return value of integral()

  // Contribution of atom j to the Born integral of atom i, as in GBn2.
  function integral(i, j, r, invR) {
    const or = born[i], sr = scaled[j];
    let value = 0;
    slope = 0;
    if (r + sr - or >= 0) {
      const d = r - sr, D = d < 0 ? -d : d, outside = D > or, L = outside ? D : or, U = r + sr;
      const iL = 1 / L, iU = 1 / U, iL2 = iL * iL, iU2 = iU * iU;
      const s2 = sr * sr * invR, w = r - s2, log = Math.log(L * iU);
      value = 0.5 * (iL - iU + 0.25 * w * (iU2 - iL2) + 0.5 * log * invR);
      slope = 0.5 * (iU2 - 0.5 * w * iU2 * iU - 0.5 * iU * invR
        + 0.25 * (1 + s2 * invR) * (iU2 - iL2) - 0.5 * log * invR * invR);
      if (outside) {
        const dL = 0.5 * (-iL2 + 0.5 * w * iL2 * iL + 0.5 * iL * invR);
        slope += d < 0 ? -dL : dL;
      }
    }
    if (r <= radius[i] + radius[j] + NECK_CUT) {
      const t = index[i] + size * index[j], z = r - d0[t], z2 = z * z, z4 = z2 * z2;
      const den = 1 / (1 + 100 * z2 + 300000 * z4 * z2), top = NECK_SCALE * m0[t] * den;
      value += top;
      slope -= top * den * (200 * z + 1800000 * z4 * z);
    }
    return value;
  }

  return function evaluate(x, g = new Float64Array(3 * n)) {
    g.fill(0);
    const E = { bond: 0, angle: 0, torsion: 0, lj_rep: 0, lj_disp: 0, coulomb: 0, gb_solv: 0 };
    bondedForces(sys, x, g, E);

    // Pass 1: Lennard-Jones, Coulomb and Born integrals.
    const first = fixedI === null;
    if (first) fixedI = new Float64Array(n);
    I.fill(0);
    skip.fill(-1);
    let rep = 0, disp = 0, coul = 0, cursor = 0;
    for (let i = 0; i < n; i++) {
      for (const j of excluded[i]) skip[j] = i;
      const a = 3 * i, xi = x[a], yi = x[a + 1], zi = x[a + 2], mi = moving[i];
      for (let j = i + 1; j < n; j++) {
        const live = mi | moving[j];
        if (!live && !first) continue;
        const b = 3 * j, dx = xi - x[b], dy = yi - x[b + 1], dz = zi - x[b + 2];
        const r2 = Math.max(dx * dx + dy * dy + dz * dz, 1e-8), r = Math.sqrt(r2), invR = 1 / r;
        const Ii = integral(i, j, r, invR), si = slope;
        const Ij = integral(j, i, r, invR), sj = slope;
        let eRep = 0, eDisp = 0, eCoul = 0;
        if (skip[j] !== i) {
          const sr = (halfSigma[i] + halfSigma[j]) * invR, sr2 = sr * sr, sr6 = sr2 * sr2 * sr2;
          const e4 = 4 * rootEps[i] * rootEps[j];
          eRep = e4 * sr6 * sr6; eDisp = -e4 * sr6; eCoul = COULOMB * q[i] * q[j] * invR;
        }
        if (live) {
          I[i] += Ii; I[j] += Ij;
          dI[cursor++] = si * invR; dI[cursor++] = sj * invR;
          rep += eRep; disp += eDisp; coul += eCoul;
          const f = -(12 * eRep + 6 * eDisp + eCoul) * invR * invR;
          g[a] += f * dx; g[a + 1] += f * dy; g[a + 2] += f * dz;
          g[b] -= f * dx; g[b + 1] -= f * dy; g[b + 2] -= f * dz;
        } else {
          fixedI[i] += Ii; fixedI[j] += Ij;
          fixedRep += eRep; fixedDisp += eDisp; fixedCoul += eCoul;
        }
      }
    }
    for (const [i, j, qq, s, e] of pairs14) {
      const a = 3 * i, b = 3 * j, dx = x[a] - x[b], dy = x[a + 1] - x[b + 1], dz = x[a + 2] - x[b + 2];
      const invR = 1 / Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), 1e-4), sr6 = (s * invR) ** 6;
      const eRep = 4 * e * sr6 * sr6, eDisp = -4 * e * sr6, eCoul = COULOMB * qq * invR;
      rep += eRep; disp += eDisp; coul += eCoul;
      const f = -(12 * eRep + 6 * eDisp + eCoul) * invR * invR;
      g[a] += f * dx; g[a + 1] += f * dy; g[a + 2] += f * dz;
      g[b] -= f * dx; g[b + 1] -= f * dy; g[b + 2] -= f * dz;
    }
    E.lj_rep = rep + fixedRep; E.lj_disp = disp + fixedDisp; E.coulomb = coul + fixedCoul;

    // Pass 2: Born radii, polar and nonpolar solvation.
    let solv = 0;
    for (let i = 0; i < n; i++) {
      const or = born[i], R = radius[i], psi = (I[i] + fixedI[i]) * or;
      const t = Math.tanh(psi * (alpha[i] + psi * (-beta[i] + psi * gamma[i])));
      const bi = 1 / (1 / or - t / (or + OFFSET));
      B[i] = bi;
      dBdI[i] = bi * bi * (1 - t * t) * (alpha[i] + psi * (-2 * beta[i] + psi * 3 * gamma[i])) * or / R;
      const self = -0.5 * GB_K * q[i] * q[i] / bi;
      const rb = R / bi, rb2 = rb * rb, ace = 28.3919551 * (R + 0.14) * (R + 0.14) * rb2 * rb2 * rb2;
      solv += self + ace;
      dEdB[i] = -(self + 6 * ace) / bi;
    }
    for (let i = 0; i < n; i++) {
      const a = 3 * i, xi = x[a], yi = x[a + 1], zi = x[a + 2], kq = GB_K * q[i], bi = B[i], mi = moving[i];
      let gx = 0, gy = 0, gz = 0, db = 0;
      for (let j = i + 1; j < n; j++) {
        const b = 3 * j, dx = xi - x[b], dy = yi - x[b + 1], dz = zi - x[b + 2];
        const r2 = dx * dx + dy * dy + dz * dz, bb = bi * B[j], u = r2 / (4 * bb);
        const ex = u > 40 ? 0 : Math.exp(-u); // below 1e-17 of f^2
        const inv = 1 / Math.sqrt(r2 + bb * ex), e = kq * q[j] * inv, e3 = e * inv * inv;
        solv -= e;
        if (ex) {
          const c = 0.5 * e3 * ex * (1 + u);
          db += c * B[j];
          dEdB[j] += c * bi;
        }
        if (mi | moving[j]) {
          const f = e3 * (1 - 0.25 * ex), fx = f * dx, fy = f * dy, fz = f * dz;
          gx += fx; gy += fy; gz += fz;
          g[b] -= fx; g[b + 1] -= fy; g[b + 2] -= fz;
        }
      }
      dEdB[i] += db;
      g[a] += gx; g[a + 1] += gy; g[a + 2] += gz;
    }
    E.gb_solv = solv;

    // Pass 3: chain rule through the Born radii.
    cursor = 0;
    for (let i = 0; i < n; i++) dEdB[i] *= dBdI[i];
    for (let i = 0; i < n; i++) {
      const a = 3 * i, xi = x[a], yi = x[a + 1], zi = x[a + 2], mi = moving[i], ei = dEdB[i];
      for (let j = i + 1; j < n; j++) {
        if (!(mi | moving[j])) continue;
        const f = ei * dI[cursor] + dEdB[j] * dI[cursor + 1];
        cursor += 2;
        const b = 3 * j, fx = f * (xi - x[b]), fy = f * (yi - x[b + 1]), fz = f * (zi - x[b + 2]);
        g[a] += fx; g[a + 1] += fy; g[a + 2] += fz;
        g[b] -= fx; g[b + 1] -= fy; g[b + 2] -= fz;
      }
    }

    E.total = E.bond + E.angle + E.torsion + E.lj_rep + E.lj_disp + E.coulomb + E.gb_solv;
    return E;
  };
}
