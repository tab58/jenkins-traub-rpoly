'use strict';

/**
 * Computes synthetic division of polynomials. Polynomial degrees n >= m.
 * @param {number[]} P Polynomial such that p[0]*x^n + p[1]*x^{n-1} + ... + p[n]*x + p[n+1]
 * @param {number[]} D Polynomial such that d[0]*x^m + d[1]*x^{m-1} + ... + d[m]*x + d[m+1]
 * @returns {Object} Computed quotient and remainder.
 */
const synthDiv = function (P, D) {
  // do input validation check
  const n = P.length - 1;
  const m = D.length - 1;
  if (m > n) {
    throw new Error('Divisor polynomial must be lower degree than dividend polynomial.');
  }
  const Q = P.slice();
  const q = n - m + 1;
  const h = D[0]; // scale for non-monic polynomials
  for (let i = 0; i < q; ++i) {
    Q[i] /= h;
    const p = Q[i];
    for (let j = 1; j <= m; ++j) {
      Q[i + j] += p * -D[j];
    }
  }
  return {
    q: Q.slice(0, q),
    r: Q.slice(q)
  };
};

module.exports = synthDiv;
