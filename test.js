'use strict';
/* globals describe it */
const expect = require('chai').expect;
const jenkinsTraub = require('./index.js');
const Complex = require('./complex.js');

const complexNumberSort = (a, b) => {
  if (a.real < b.real ||
      a.imag < b.imag) {
    return -1;
  }
  if (a.real > b.real ||
      a.imag > b.imag) {
    return 1;
  }
  return 0;
};

describe('Jenkins-Traub Roots', () => {
  it('Test 1 (from paper)', () => {
    // TODO: Compute real bound for roots, not just this magic number
    const TOL = 1e-7;

    const P = [1, -6.01, 12.54, -8.545, -5.505, 12.545, -8.035, 2.01];
    const roots = jenkinsTraub(P).sort(complexNumberSort);
    const verifiedRoots = [
      new Complex(0.5, 0.5),
      new Complex(0.5, -0.5),
      new Complex(1, 0),
      new Complex(1, 0),
      new Complex(-1, 0),
      new Complex(2, 0),
      new Complex(2.01, 0)
    ].sort(complexNumberSort);

    expect(roots.length).to.be.eql(verifiedRoots.length);
    for (let i = 0; i < roots.length; ++i) {
      const root = roots[i];
      const verifiedRoot = verifiedRoots[i];
      expect(Math.abs(root.real - verifiedRoot.real)).to.be.below(TOL);
      expect(Math.abs(root.imag - verifiedRoot.imag)).to.be.below(TOL);
      expect(Math.abs(root.abs() - verifiedRoot.abs())).to.be.below(TOL);
    }
  });
});
