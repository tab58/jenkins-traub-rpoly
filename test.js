'use strict';
/* globals describe it */
const expect = require('chai').expect;
const jenkinsTraub = require('./index.js');

describe('Jenkins-Traub Roots', () => {
  it('Test 1 (from paper)', () => {
    const P = [1, -6.01, 12.54, -8.545, -5.505, 12.545, -8.035, 2.01];
    const roots = jenkinsTraub(P);
    const realRoots = [
      { real: 0.5, imag: 0.5 },
      { real: 0.5, imag: -0.5 },
      { real: 1, imag: 0 },
      { real: 1, imag: 0 },
      { real: -1, imag: 0 },
      { real: 2, imag: 0 },
      { real: 2.01, imag: 0 }
    ];
    expect(roots).to.be.eql(realRoots);
  });
});
