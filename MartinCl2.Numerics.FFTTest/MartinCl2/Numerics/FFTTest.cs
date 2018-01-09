using System;
using System.Collections.Generic;
using System.Numerics;
using Xunit;

namespace MartinCl2.Numerics
{
    public class FFTTest
    {
        private static double PRECISION = 1e-7;

        private class ComplexComparer : IEqualityComparer<Complex>
        {
            public bool Equals(Complex x, Complex y)
            {
                return (x - y).Magnitude < PRECISION;
            }

            public int GetHashCode(Complex obj)
            {
                throw new NotImplementedException();
            }
        }

        protected static IEqualityComparer<Complex> comparer = new ComplexComparer();

        [Fact]
        public void SinglePointTest()
        {
            FFTCalculator fft = new FFTCalculator(1);
            Complex[] signal = new Complex[1];
            Complex[] spectrum = new Complex[1];
            Complex[] ifft = new Complex[1];
            signal[0] = new Complex(1, 1);
            fft.DFT(signal, spectrum);
            Assert.Equal(signal, spectrum, comparer);
            fft.IDFT(spectrum, ifft);
            Assert.Equal(signal, ifft, comparer);
        }

        [Fact]
        public void DualPointTest()
        {
            FFTCalculator fft = new FFTCalculator(2);

            Complex[] signal = new Complex[2];
            signal[0] = new Complex(1.7, 3.2);
            signal[1] = new Complex(-2.3, 1.6);

            Complex[] expected = new Complex[2];
            expected[0] = new Complex(-0.6, 4.8);
            expected[1] = new Complex(4, 1.6);

            Complex[] spectrum = new Complex[2];
            fft.DFT(signal, spectrum);
            Assert.Equal(expected, spectrum, comparer);

            Complex[] ifft = new Complex[2];
            fft.IDFT(spectrum, ifft);
            Assert.Equal(signal, ifft, comparer);
        }

        [Fact]
        public void QuadPointTest()
        {
            FFTCalculator fft = new FFTCalculator(4);

            Complex[] signal = new Complex[4];
            signal[0] = new Complex(0.3, -2.7);
            signal[1] = new Complex(-2.2, -1.8);
            signal[2] = new Complex(3.1, -2.5);
            signal[3] = new Complex(3.1, 0.7);

            Complex[] expected = new Complex[4];
            expected[0] = new Complex(4.3, -6.3);
            expected[1] = new Complex(-5.3, 5.1);
            expected[2] = new Complex(2.5, -4.1);
            expected[3] = new Complex(-0.3, -5.5);

            Complex[] spectrum = new Complex[4];
            fft.DFT(signal, spectrum);
            Assert.Equal(expected, spectrum, comparer);

            Complex[] ifft = new Complex[4];
            fft.IDFT(spectrum, ifft);
            Assert.Equal(signal, ifft, comparer);
        }
    }
}
