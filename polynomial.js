'use strict';

const Complex = require('./complex.js');

/**
 * Computes the derivative of a polynomial.
 * @param {number[]} P The polynomial to compute the derivative.
 * @returns {number[]} The derivative of the polynomial P.
 */
const computePolyDeriv = function (P) {
  const nn = P.length - 1;
  const D = [];
  for (let i = 0; i < nn; ++i) {
    D.push((nn - i) * P[i]);
  }
  return D;
};

/**
 * Evaluates a polynomial at the given point.
 * @param {number[]} P The polynomial to evaluate.
 * @param {number} x The evaluation point of the polynomial.
 * @returns {number} The value of the polynomial at the evaluation point.
 */
const hornerPolyEval = function (P, x) {
  const nn = P.length - 1;
  let ff = P[0];
  for (let i = 1; i <= nn; ++i) {
    ff = ff * x + P[i];
  }
  return ff;
};

/**
 * Evaluates the polynomial and its derivative at the given point.
 * @param {number[]} P The polynomial.
 * @param {number} x The value at which the polynomial will be evaluated.
 * @param {Object} out An object that will contains the output values (optional).
 * @returns {Object} The out argument is specified, a new object if not.
 */
const hornerPolyDerivEval = function (P, x, out) {
  const nn = P.length - 1;
  const n = nn - 1;
  let ff = P[0];
  let df = ff;
  for (let i = 1; i <= n; ++i) {
    ff = ff * x + P[i];
    df = df * x + ff;
  }
  ff = ff * x + P[nn];

  const V = out || { f: 0, df: 0 };
  V.f = ff;
  V.df = df;
  return V;
};

/**
 * Evaluates a complex polynomial with all real coefficients.
 * @param {number[]} P The polynomial.
 * @param {Complex} cX The complex evaluation point.
 * @returns {Complex} The value of the polynomial at the evaluation point.
 */
const hornerComplexPolyEval = function (P, cX) {
  const x = cX.real;
  const y = cX.imag;
  let real = P[0];
  let imag = 0;

  for (let i = 1; i < P.length; ++i) {
    const a = P[i];
    const u = real;
    const v = imag;
    real = u * x - v * y + a;
    imag = v * x + u * y;
  }
  return new Complex(real, imag);
};

const scalePolynomial = function (P, scale) {
  const n = P.length - 1;
  for (let i = 0; i <= n; ++i) {
    P[i] *= scale;
  }
  return P;
};

const copyPolynomial = function (copyFrom, copyTo) {
  const n = copyFrom.length;
  if (n !== copyTo.length) {
    throw new Error('copyPolynomial(): Polynomial arrays must be same length.');
  }
  for (let i = 0; i < n; ++i) {
    copyTo[i] = copyFrom[i];
  }
};

/**
 * Adds an array of polynomials.
 * @param {number[][]} polyArray An array of polynomials.
 * @returns {number[]} The sum of the polynomials.
 */
const addPolynomials = function (polyArr, polyScales) {
  if (polyArr.length !== polyScales.length) {
    throw new Error('addPolynomials(): scale array must be the same length as the polynomial array.');
  }
  const initLen = polyArr.length > 0 ? polyArr[0].length : 0;
  if (initLen === 0) {
    return undefined;
  }
  const idx = [];
  for (let j = 0; j < polyArr.length; ++j) {
    idx.push(j);
  }
  const sortedIndices = idx.sort((a, b) => polyArr[b].length - polyArr[a].length);
  const alen = sortedIndices.length;
  const firstIdx = sortedIndices[0];
  const S = polyArr[firstIdx].slice();
  const scl = polyScales ? polyScales[firstIdx] : 1;
  for (let j = 0; j < S.length; ++j) {
    S[j] *= scl;
  }

  const sn = S.length - 1;
  for (let j = 1; j < alen; ++j) {
    const idx = sortedIndices[j];
    const P = polyArr[idx];
    const scale = polyScales ? polyScales[idx] : 1;
    const pn = P.length - 1;
    for (let i = 0; i <= pn; ++i) {
      S[sn - i] += scale * P[pn - i];
    }
  }
  return S;
};

module.exports = {
  hornerComplexPolyEval,
  hornerPolyDerivEval,
  hornerPolyEval,
  computePolyDeriv,
  scalePolynomial,
  addPolynomials,
  copyPolynomial
};
