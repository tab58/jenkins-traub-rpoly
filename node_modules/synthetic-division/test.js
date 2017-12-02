'use strict';
/* globals describe it */
const expect = require('chai').expect;

const synthDiv = require('./index.js');

describe('Synthetic Division', () => {
  it('Degree 1 Divisor', () => {
    const P = [1, -12, 0, -42];
    const D = [1, -3];
    const { q: Q, r: R } = synthDiv(P, D);
    expect(Q).to.be.eql([1, -9, -27]);
    expect(R).to.be.eql([-123]);
  });
  it('Degree 2 Divisor', () => {
    const P = [1, -12, 0, -42];
    const D = [1, 1, -3];
    const { q: Q, r: R } = synthDiv(P, D);
    expect(Q).to.be.eql([1, -13]);
    expect(R).to.be.eql([16, -81]);
  });
  it('Non-monic Polynomials', () => {
    const P = [6, 5, 0, -7];
    const D = [3, -2, -1];
    const { q: Q, r: R } = synthDiv(P, D);
    expect(Q).to.be.eql([2, 3]);
    expect(R).to.be.eql([8, -4]);
  });
});
