'use strict';
/* globals describe it */
const expect = require('chai').expect;
const jenkinsTraub = require('./index.js');
const Complex = require('./complex.js');

// sorts the complex numbers by ascending reals, then ascending imaginary components.
const complexNumberSort = (a, b) => {
  if (a.real < b.real) {
    return -1;
  }
  if (a.real > b.real) {
    return 1;
  }
  if (a.imag < b.imag) {
    return -1;
  }
  if (a.imag > b.imag) {
    return 1;
  }
  return 0;
};

const expectToBeNearlyEqual = (a, b, TOL) => {
  return expect(Math.abs(a - b)).to.be.below(TOL);
};

const checkPolynomialRoots = (P, roots, TOL) => {
  const calculatedRoots = jenkinsTraub(P).sort(complexNumberSort);
  const verifiedRoots = roots.sort(complexNumberSort);
  expect(calculatedRoots.length).to.be.eql(verifiedRoots.length);
  for (let i = 0; i < calculatedRoots.length; ++i) {
    const root = calculatedRoots[i];
    const verifiedRoot = verifiedRoots[i];
    const TOLERANCE = Array.isArray(TOL) ? TOL[i] : TOL;
    expectToBeNearlyEqual(root.real, verifiedRoot.real, TOLERANCE);
    expectToBeNearlyEqual(root.imag, verifiedRoot.imag, TOLERANCE);
    expectToBeNearlyEqual(root.abs(), verifiedRoot.abs(), TOLERANCE);
  }
};

// TODO: Compute real bounds for roots, not just these magic numbers
describe('Jenkins-Traub Roots', () => {
  it('Smoke Test - from theoretical paper', () => {
    const P = [1, -6.01, 12.54, -8.545, -5.505, 12.545, -8.035, 2.01];
    const verifiedRoots = [
      new Complex(0.5, 0.5),
      new Complex(0.5, -0.5),
      new Complex(1, 0),
      new Complex(1, 0),
      new Complex(-1, 0),
      new Complex(2, 0),
      new Complex(2.01, 0)
    ];
    checkPolynomialRoots(P, verifiedRoots, 1e-7);
  });
  describe('Test Cubics', () => {
    describe('Small Integer Coefficients', () => {
      const GEN_TOL = 1e-10;
      it('Test 1', () => {
        checkPolynomialRoots([1, -6, 11, -6], [
          new Complex(3, 0),
          new Complex(1, 0),
          new Complex(2, 0)
        ], GEN_TOL);
      });
      it('Test 2', () => {
        checkPolynomialRoots([1, 0, 0, 1], [
          new Complex(-1, 0),
          new Complex(0.5, Math.sqrt(0.75)),
          new Complex(0.5, -Math.sqrt(0.75))
        ], GEN_TOL);
      });
      it('Test 3', () => {
        checkPolynomialRoots([1, -3, 2, 0], [
          new Complex(0, 0),
          new Complex(1, 0),
          new Complex(2, 0)
        ], GEN_TOL);
      });
      it('Test 4', () => {
        checkPolynomialRoots([1, -30, 299, -1980], [
          new Complex(20, 0),
          new Complex(5, Math.sqrt(74)),
          new Complex(5, -Math.sqrt(74))
        ], GEN_TOL);
      });
    });

    describe('Zeros of Different Magnitudes', () => {
      const GEN_TOL = 1e-10;
      const t = 1e-5;
      const h = 1e9;
      it('Test 1', () => {
        checkPolynomialRoots([1, -30, 299, -t], [
          new Complex(t / 299, 0),
          new Complex(15, Math.sqrt(74)),
          new Complex(15, -Math.sqrt(74))
        ], [
          GEN_TOL,
          1e-7,
          1e-7
        ]);
      });
      it('Test 2', () => {
        checkPolynomialRoots([t, -h, h, -t], [
          new Complex(1, 0),
          new Complex(t / h, 0),
          new Complex(h / t, 0)
        ], [
          GEN_TOL,
          GEN_TOL,
          1
        ]);
      });
      it('Test 3', () => {
        checkPolynomialRoots([1, -h, -t, h * t], [
          new Complex(h, 0),
          new Complex(Math.sqrt(t), 0),
          new Complex(-Math.sqrt(t), 0)
        ], [
          1,
          GEN_TOL,
          GEN_TOL
        ]);
      });
    });

    describe('Ill-conditioned Zeros', () => {
      const GEN_TOL = 1e-10;
      const t = 1e-5;
      const h = 1e9;
      const M = Math.pow(2, 5);
      const N = Math.pow(2, 40) - 1;
      const u = M / N;
      const v = 1 / (2 * N);
      it('Test 1', () => {
        checkPolynomialRoots([N + 1, -(N - 1), -(N + 1), N - 1], [
          new Complex(-1, 0),
          new Complex(1, 0),
          new Complex(1 - 2 / (N + 1), 0)
        ], [
          GEN_TOL,
          GEN_TOL,
          1e-2
        ]);
      });
      // it('Test 2', () => {
      //   const C = 9 * Math.pow(N, 3);
      //   checkPolynomialRoots([3 * C, 3 * C, C, 1 - Math.pow(N, 3)], [
      //     new Complex((1 - 2 * v) / 3, 0),
      //     new Complex(1 + v, v * Math.sqrt(3)),
      //     new Complex(1 + v, -v * Math.sqrt(3))
      //   ], [
      //     GEN_TOL,
      //     GEN_TOL,
      //     1e-2
      //   ]);
      // });
    });
  });
});
