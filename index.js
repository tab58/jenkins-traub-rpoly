'use strict';

const _Math = Math;
const {
  hornerPolyEval,
  hornerComplexPolyEval,
  scalePolynomial,
  computePolyDeriv,
  copyPolynomial,
  addPolynomials
} = require('./polynomial.js');
const Complex = require('./complex.js');
const synthDiv = require('synthetic-division');
const qdrtc = require('./quadRoots.js');

const DEG2RAD = _Math.PI / 180;
const XX = _Math.cos(94 * DEG2RAD);
const YY = _Math.sin(94 * DEG2RAD);
const X0 = _Math.cos(49 * DEG2RAD);
const Y0 = _Math.sin(49 * DEG2RAD);

/**
 * Computes a single Newton-Raphson iteration from the guess.
 * @param {number[]} P The real polynomial.
 * @param {number} x The initial guess.
 * @returns {number} The new guess calculated via Newton-Raphson.
 */
const newtonsPolyIteration = function (P, x) {
  const nn = P.length - 1;
  const n = nn - 1;
  let ff = P[0];
  let df = ff;
  for (let i = 1; i <= n; ++i) {
    ff = ff * x + P[i];
    df = df * x + ff;
  }
  ff = ff * x + P[nn];
  return x - (ff / df);
};

/**
 * Computes the Cauchy lower bound of the set of roots of the given polynomial.
 * @param {number[]} P The given polynomial.
 */
const getRootBound = function (P) {
  const nn = P.length - 1;
  const n = nn - 1;
  const pt = P.map(val => _Math.abs(val));
  pt[nn] = -pt[nn];
  // get a guess
  let x = _Math.pow(_Math.E, (_Math.log(-pt[nn]) - _Math.log(pt[0])) / n);
  // check the coefficients to get the lowest estimate
  if (pt[n] !== 0) {
    const xm = -pt[nn] / pt[n];
    if (xm < x) {
      x = xm;
    }
  }
  // cut the interval (0,x) down by half
  let ff = 0;
  do {
    x /= 2;
    ff = hornerPolyEval(pt, x);
  } while (ff > 0);
  // do Newton-Raphson until we get 2 decimal places
  let x1 = x;
  do {
    x = x1;
    x1 = newtonsPolyIteration(pt, x);
  } while (_Math.abs(1 - (x1 / x)) > 0.005);
  return x1;
};

/**
 * Computes a no-shift sequence of polynomials. Requires the first coefficient of K0 to be 1.
 * @param {number[]} P The original polynomial in the sequence.
 * @param {number} M The number of the polynomial in the sequence to evaluate.
 * @returns {number[]} The Mth polynomial in the sequence.
 */
const computeNoShiftK = function (K0, P, M) {
  const pn = P.length - 1;
  const n = K0.length - 1;
  const P0 = P[pn]; // P(0)
  if (P0 === 0) {
    throw new Error('Constant term in P should not be zero.');
  }
  // do no-shift
  const K = K0.slice();
  for (let j = 0; j < M; ++j) {
    const Klam0 = K[n]; // K(lam)(0)
    const t = -Klam0 / P0;
    for (let i = n; i >= 1; --i) {
      K[i] = K[i - 1] + t * P[i];
    }
    K[0] = t * P[0]; // assumes P[0] = 1
  }
  return K;
};

/**
 * Computes the next K^(lambda+1)(z).
 * @param {number[]} K0 The previous K^(lambda) polynomial in the sequence.
 * @param {number[]} sigma The quadratic factor used in the fixed shift.
 * @param {number[]} Qp The quotient of P(z) / sigma(z).
 * @param {number[]} Qk The quotient of K^(lambda)(z) / sigma(z).
 * @param {number} a P(s1) = a - b*s2
 * @param {number} b P(s1) = a - b*s2
 * @param {number} c K^(lambda)(s1) = c - d*s2
 * @param {number} d K^(lambda)(s1) = a - d*s2
 */
const computeNextFixedShiftK = function (K0, Qp, Qk, a, b, c, d, u, v) {
  const K = K0.slice();
  const n = K.length - 1;
  const nqp = Qp.length - 1;
  const nqk = Qk.length - 1;
  // compute new K^(lambda+1) parameters
  const alpha = a * a + u * a * b + v * b * b;
  const beta = -(a * c + u * a * d + v * b * d);
  const gamma = b * c - a * d;
  // determines whether to use the scaled or regular recurrence
  const useScaled = (_Math.abs(gamma) > 1e-15);
  const Qkc = useScaled ? (alpha / gamma) : 1;
  const Qzc = useScaled ? 1 : (gamma / alpha);
  const Qpc = useScaled ? (beta / gamma) : (beta / alpha);
  for (let i = 0; i <= n; ++i) {
    let ki = 0;
    if (nqk - i >= 0) {
      ki += Qkc * Qk[nqk - i];
    }
    if (nqp - i >= 0) {
      ki += Qpc * Qp[nqp - i];
    }
    if (i - 1 >= 0) {
      ki += Qzc * Qp[nqp - i + 1];
    }
    K[n - i] = ki;
  }
  K[n] += Qzc * b;
  return K;
};

/**
 * Computes the next approximate sigma with parameters. Method detailed in "Three Stage
 * Variable-Shift Iterations for the Solution of Polynomial Equations With a Posteriori
 * Error Bounds for the Zeros" by M.A. Jenkins, Doctoral Thesis, Stanford University, 1969.
 * @param {number} a P(s1) = a - b * s2.
 * @param {number} b P(s1) = a - b * s2.
 * @param {number} c K^{lambda}(s1) = c - d * s2.
 * @param {number} d K^{lambda}(s1) = c - d * s2.
 * @param {number} u sigma(z) = z^z + u * z + v.
 * @param {number} v sigma(z) = z^z + u * z + v.
 * @param {number} alpha K^{lambda+1}(z) = 1/z * (K^{lambda}(z) - alpha * P(z)).
 */
const computeSigmaEstimate = function (a, b, c, d, u, v, K, P) {
  const kn = K.length - 1;
  const pn = P.length - 1;
  const alpha0 = -K[kn] / P[pn];
  const alpha1 = -(K[kn - 1] + alpha0 * P[pn - 1]) / P[pn];

  const a1 = b * c - a * d;
  const a2 = a * c + u * a * d + v * b * d;
  const c2 = alpha0 * a2;
  const c3 = alpha0 * alpha0 * (a * a + u * a * b + v * b * b);
  const c4 = v * alpha1 * a1 - c2 - c3;
  const c1 = c * c + u * c * d + v * d * d + alpha0 * (a * c + u * b * c + v * b * d) - c4;
  const dUNum = -(u * (c2 + c3) + v * (alpha0 * a1 + alpha1 * a2));
  const dVNum = v * c4;

  const deltaU = _Math.abs(dUNum) < 1e-15 && _Math.abs(c1) < 1e-15 ? 0 : dUNum / c1;
  const deltaV = _Math.abs(dVNum) < 1e-15 && _Math.abs(c1) < 1e-15 ? 0 : dVNum / c1;

  // Update u and v in the quadratic sigma.
  return [1, u + deltaU, v + deltaV];
};

/**
 * Determines if the factor converges to a linear factor.
 * @param {number[]} arr Array of linear root approximations.
 * @param {*} i The last index of the array (index is i mod 3).
 */
const hasLinearConverged = function (arr, i) {
  if (i >= 2) {
    const inext = i % 3;
    const icurr = (i - 1) % 3;
    const iprev = (i - 2) % 3;
    const tnext = arr[inext];
    const tcurr = arr[icurr];
    const tprev = arr[iprev];
    return tnext.clone().sub(tcurr).abs() <= 0.5 * tcurr.abs() &&
           tcurr.clone().sub(tprev).abs() <= 0.5 * tprev.abs();
  }
  return false;
};

/**
 * Determines if the factor converges to a quadratic factor.
 * @param {number[]} arr Array of quadratic root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasQuadraticConverged = function (arr, i) {
  if (i >= 2) {
    const inext = i % 3;
    const icurr = (i - 1) % 3;
    const iprev = (i - 2) % 3;
    const vnext = arr[inext];
    const vcurr = arr[icurr];
    const vprev = arr[iprev];
    return _Math.abs(vnext - vcurr) <= 0.5 * _Math.abs(vcurr) &&
           _Math.abs(vcurr - vprev) <= 0.5 * _Math.abs(vprev);
  }
  return false;
};

/**
 * Establishes criteria to tell if a real root has converged.
 * @param {number[]} roots The array of root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasRealRootConverged = function (roots, i) {
  if (i >= 2) {
    const r2 = roots[i % 3];
    const r1 = roots[(i - 1) % 3];
    const r0 = roots[(i - 2) % 3];
    const ei = _Math.abs(r2 - r1);
    const ei1 = _Math.abs(r1 - r0);
    const magRoot = _Math.abs(r1);
    return wardCriterion(ei, ei1, magRoot);
  }
  return false;
};

/**
 * Establishes criteria to tell if a complex root has converged.
 * @param {number[]} roots The array of root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasComplexRootConverged = function (croots, i) {
  if (i >= 2) {
    const r2 = croots[i % 3];
    const r1 = croots[(i - 1) % 3];
    const r0 = croots[(i - 2) % 3];
    const ei = r2.clone().sub(r1).abs();
    const ei1 = r1.clone().sub(r0).abs();
    const magRoot = r1.abs();
    return wardCriterion(ei, ei1, magRoot);
  }
  return false;
};

// /**
//  * The absolute root tolerance, i.e. if P(z) is below this, it's a root.
//  */
// const ABSOLUTE_ROOT_TOLERANCE = 1e-15;

// /**
//  * The relative root tolerance. Used to test for convergence.
//  */
// const RELATIVE_ROOT_TOLERANCE = 1e-10;

// /**
//  * The tolerance of magnitude of the root for it to be considered a root.
//  */
// const ROOT_MAGNITUDE_TOLERANCE = 1e-4;

/**
 * Ward's stopping criterion as described in "New stopping criteria for iterative root finding"
 * by Nikolajsen, Jorgen L., Royal Society open science (2014).
 * @param {*} ei The absolute magnitude of root separation between the most recent approximations, i.e. |z_i - z_(i-1)|.
 * @param {*} ei1 The absolute magnitude of root separation between the older calculations, i.e. |z_(i-1) - z_(i-2)|.
 * @param {*} magRoot The magnitude of the next to last root, i.e. |z_(i-1)|.
 */
const wardCriterion = function (ei, ei1, magRoot) {
  if ((magRoot < 1e-4 && ei <= 1e-7) ||
      (magRoot >= 1e-4 && ei / magRoot <= 1e-3)) {
    return ei >= ei1;
  }
  return false;
  // if (ei <= ei1) {
  //   if (magRoot < ROOT_MAGNITUDE_TOLERANCE) {
  //     return ei < ABSOLUTE_ROOT_TOLERANCE;
  //   } else {
  //     return ei / magRoot <= RELATIVE_ROOT_TOLERANCE;
  //   }
  // }
  // return false;
};

/**
 * Computes Stage 3 of the Jenkins-Traub algorithm for convergence to a linear factor.
 * @param {number[]} K0 The last K polynomial calculated in Stage 2.
 * @param {number[]} P The polynomial for which the roots are to be solved.
 * @param {Complex} sL The last root approximation calculated in Stage 2.
 * @param {number} L The number of Stage 2 iterations.
 */
const computeLinearVariableShift = function (K0, P, sL, L) {
  const n = K0.length - 1;
  let sReal = sL.clone().sub(hornerComplexPolyEval(P, sL).divide(hornerComplexPolyEval(K0, sL))).real;
  let K = K0.slice();
  const factor = [1, -sReal];
  const rootApproximations = [0, 0, 0];
  for (let k = 0; k < 10 * L; ++k) {
    factor[1] = -sReal;
    rootApproximations[k % 3] = sReal;
    // construct next K similar to no-shift sequence
    const { q: Qp, r: rP } = synthDiv(P, factor);
    const { q: Qk, r: rK } = synthDiv(K, factor);
    let KsReal = rK[0];
    let PsReal = rP[0];
    let t = -KsReal / PsReal;
    for (let i = n; i >= 1; --i) {
      K[i] = Qk[i - 1] + t * Qp[i];
    }
    K[0] = t * Qp[0]; // assumes P[0] = 1
    sReal -= PsReal * K[0] / hornerPolyEval(K, sReal);

    // TODO: find convergence criteria
    if (hasRealRootConverged(rootApproximations, k)) {
      const root = rootApproximations[(k - 1) % 3];
      return {
        roots: [new Complex(root, 0)]
      };
    }
  }
  return undefined;
};

/**
 * Computes Stage 3 of the Jenkins-Traub algorithm for convergence to a quadratic factor.
 * @param {number[]} K0 The last K polynomial calculated in Stage 2.
 * @param {number[]} P The polynomial for which the roots are to be solved.
 * @param {Complex} sL The last root approximation calculated in Stage 2.
 * @param {number} L The number of Stage 2 iterations.
 */
const computeQuadraticVariableShift = function (K0, P, sigmaLambda, L) {
  let K = K0.slice();
  let sigmaL = sigmaLambda.slice();
  let uLambda = sigmaL[1];
  let vLambda = sigmaL[2];
  const rootApproximations = [];
  const s = new Complex(0, 0);
  for (let i = 0; i < 10 * L; ++i) {
    const { q: Qp, r: remP } = synthDiv(P, sigmaL);
    const b = remP.length > 1 ? remP[0] : 0;
    const a = (remP.length > 1 ? remP[1] : remP[0]) - b * uLambda;
    const { q: Qk, r: remK } = synthDiv(K, sigmaL);
    const d = remK.length > 1 ? remK[0] : 0;
    const c = (remK.length > 1 ? remK[1] : remP[0]) - d * uLambda;

    const [, uLambda1, vLambda1] = computeSigmaEstimate(a, b, c, d, uLambda, vLambda, K, P);
    const K1 = computeNextFixedShiftK(K, Qp, Qk, a, b, c, d, uLambda1, vLambda1);
    const [sx1, sy1] = qdrtc(1, uLambda1, vLambda1);
    s.real = sx1;
    s.imag = sy1;
    rootApproximations[i % 3] = s.clone();
    if (hasComplexRootConverged(rootApproximations, i)) {
      // const Ps1 = hornerComplexPolyEval(P, rootApproximations[i % 3]).abs();
      // const Ps2 = hornerComplexPolyEval(P, rootApproximations[(i - 1) % 3]).abs();
      const root = rootApproximations[(i - 1) % 3];
      return {
        roots: [root.clone(), root.clone().conj()],
        sigma: sigmaL
      };
    }
    K = K1;
    uLambda = uLambda1;
    vLambda = vLambda1;
    sigmaL[1] = uLambda;
    sigmaL[2] = vLambda;
  }
  return undefined;
};

module.exports = function jenkinsTraub (OP) {
  if (OP[0] === 0) {
    throw new Error('Leading coefficient must not be zero.');
  }
  const zeros = [];
  // determine zeros at zero
  let n = OP.length - 1;
  while (OP[n] === 0) {
    zeros.push(new Complex(0, 0));
    n--;
  }
  // make monic polynomial
  let P = OP.slice(0, n + 1);
  scalePolynomial(P, 1 / P[0]);
  // number of roots left
  n = P.length - 1;
  const STAGE2_LIMIT = 20 * n;
  while (n > 2) {
    // Do stage 1 of algorithm for separation of zeros
    const M = 5;
    const K0 = computePolyDeriv(P);
    const KM = computeNoShiftK(K0, P, M);

    // Do stage 2-3 of the algorithm
    const rootRadius = getRootBound(P);
    let rootFound = false;
    let k = 0;
    let roots, polyFactor;
    while (!rootFound) {
      // get the base K-polynomial
      let KL = KM.slice();
      // choose a root on the radius
      const x = rootRadius * (X0 + k * XX);
      const y = rootRadius * (Y0 + k * YY);
      // const x = 0.042019;
      // const y = 0.1836611;
      const s = new Complex(x, y);
      const u = -2 * x;
      const v = x * x + y * y;
      const sigma = [1, u, v];
      const { q: Qp, r: remP } = synthDiv(P, sigma);
      const b = remP.length > 1 ? remP[0] : 0;
      const a = (remP.length > 1 ? remP[1] : remP[0]) - b * u;

      // Stage 2
      // -------
      let stage2Success = false;
      let convergingToLinear = false;
      let convergingToQuadratic = false;
      let stage2Iters = 0;
      const sigmaLambda = [1, u, v];

      while (!stage2Success) {
        const pAtRoot = new Complex(a - b * s.real, b * s.imag); // a - b * s.conj()
        let tLambdas = [];
        let vLambdas = [];
        let c, d, Qk, remK;
        // start iteration of the K-polynomials
        for (let i = 0; i <= STAGE2_LIMIT; ++i) {
          scalePolynomial(KL, 1.0 / KL[0]);
          ({ q: Qk, r: remK } = synthDiv(KL, sigma));
          d = remK.length > 1 ? remK[0] : 0;
          c = (remK.length > 1 ? remK[1] : remK[0]) - d * u;
          // get linear termination criteria
          const kAtRoot = new Complex(c - d * s.real, d * s.imag); // c - d * s2
          const t = s.clone().sub(pAtRoot.clone().divide(kAtRoot));
          tLambdas[i % 3] = t;
          convergingToLinear = hasLinearConverged(tLambdas, i);

          // get quadratic termination criteria
          const [, uLambda, vLambda] = computeSigmaEstimate(a, b, c, d, u, v, KL, P);
          vLambdas[i % 3] = vLambda;
          convergingToQuadratic = hasQuadraticConverged(vLambdas, i);

          sigmaLambda[1] = uLambda;
          sigmaLambda[2] = vLambda;
          if (convergingToLinear || convergingToQuadratic) {
            stage2Success = true;
            stage2Iters = i;
            break;
          }

          // compute new K^(lambda+1) parameters
          // remainders are either 1- or 2-length arrays
          // remP[0, 1] = [b, a + b * u]
          // remK[0, 1] = [d, c + d * u]
          if (!stage2Success) {
            KL = computeNextFixedShiftK(KL, Qp, Qk, a, b, c, d, u, v);
          }
        }
      }

      // Stage 3
      // -------
      if (convergingToQuadratic) {
        const quadInfo = computeQuadraticVariableShift(KL, P, sigmaLambda, stage2Iters);
        if (quadInfo) {
          rootFound = true;
          const { sigma, roots: quadRoots } = quadInfo;
          roots = quadRoots;
          polyFactor = sigma;
        }
      } else if (convergingToLinear) {
        const linInfo = computeLinearVariableShift(KL, P, s, stage2Iters);
        if (linInfo) {
          rootFound = true;
          const { roots: linRoot } = linInfo;
          roots = linRoot;
          polyFactor = [1, linRoot[0].real];
        }
      } else {
        throw new Error('Indeterminate state.');
      }

      if (!rootFound) {
        // increment the clocking in the complex plane and try again
        k++;
      }
    }
    // root should be found
    zeros.push.apply(zeros, roots);
    n -= roots.length;
    ({ q: P } = synthDiv(P, polyFactor));
  }
  return zeros;
};
