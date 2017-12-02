'use strict';

const _Math = Math;

/**
 * A class that represents a complex number.
 */
class Complex {
  /**
   * Creates a complex number.
   * @param {number} real The real part of the complex number.
   * @param {number} imag The imaginary part of the complex number.
   */
  constructor (real, imag) {
    this.real = real;
    this.imag = imag;
  }

  /**
   * Scale the complex number by a factor.
   * @param {number} s The scale factor.
   */
  scale (s) {
    this.real *= s;
    this.imag *= s;
    return this;
  }

  /**
   * Clones the complex number by creating another Complex object with the same values.
   */
  clone () {
    return new Complex(this.real, this.imag);
  }

  /**
   * Adds a complex number to this one.
   * @param {Complex} cX The complex number to add.
   */
  add (cX) {
    this.real += cX.real;
    this.imag += cX.imag;
    return this;
  }

  /**
   * Subtracts a complex number from this one.
   * @param {Complex} cX The complex number to subtract.
   */
  sub (cX) {
    this.real -= cX.real;
    this.imag -= cX.imag;
    return this;
  }

  /**
   * Multiplies a complex number by this one.
   * @param {Complex} cX The complex number to multiply.
   */
  multiply (cX) {
    const a = this.real;
    const b = this.imag;
    const c = cX.real;
    const d = cX.imag;

    this.real = a * c - b * d;
    this.imag = b * c + a * d;
    return this;
  }

  /**
   * Divides this complex number by a given number.
   * @param {Complex} cX The complex number to divide by.
   */
  divide (cX) {
    const a = this.real;
    const b = this.imag;
    const c = cX.real;
    const d = cX.imag;

    const den = c * c + d * d;
    this.real = (a * c + b * d) / den;
    this.imag = (b * c - a * d) / den;
    return this;
  }

  /**
   * Computes the magnitude of the complex number.
   * @returns {number} The magnitude of the complex number.
   */
  abs () {
    // deals with overflow/underflow
    const a = this.real;
    const b = this.imag;
    if (a === 0 && b === 0) {
      return 0;
    }
    const x = _Math.abs(a);
    const y = _Math.abs(b);
    const u = _Math.max(x, y);
    const t = _Math.min(x, y) / u;
    return u * _Math.sqrt(1 + t * t);
  }

  conj () {
    this.imag = -this.imag;
    return this;
  }
}

module.exports = Complex;
