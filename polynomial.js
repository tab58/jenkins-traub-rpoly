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

/**
 * Scales a polynomial in-place.
 * @param {number[]} P The polynomial to scale.
 * @param {*} scale The scalar quantity to scale the coefficients by.
 * @returns {number[]} The argument P.
 */
const scalePolynomial = function (P, scale) {
  const n = P.length - 1;
  for (let i = 0; i <= n; ++i) {
    P[i] *= scale;
  }
  return P;
};

module.exports = {
  hornerComplexPolyEval,
  hornerPolyEval,
  computePolyDeriv,
  scalePolynomial
};
