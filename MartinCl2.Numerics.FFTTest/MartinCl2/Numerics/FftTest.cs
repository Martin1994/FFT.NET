using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using Xunit;

namespace MartinCl2.Numerics
{
    public class NormalFftTest : FftTest
    {
        protected override IFftCalculator CreateCalculator(int size)
        {
            return new FftCalculator(size);
        }
    }

    public class ParallelFftTest : FftTest
    {
        protected override IFftCalculator CreateCalculator(int size)
        {
            return new ParallelFftCalculator(size);
        }
    }

    public class SpanFftTest : FftTest
    {
        protected override IFftCalculator CreateCalculator(int size)
        {
            return new SpanFftCalculator(size);
        }
    }

    public class VectorTFftTest : FftTest
    {
        protected override IFftCalculator CreateCalculator(int size)
        {
            return new VectorTFftCalculator(size);
        }
    }

    public class IntrinsicsFftTest : FftTest
    {
        protected override IFftCalculator CreateCalculator(int size)
        {
            return new IntrinsicsFftCalculator(size);
        }
    }

    public abstract class FftTest
    {
        private static double PRECISION = 1e-7;

        protected abstract IFftCalculator CreateCalculator(int size);

        private class ComplexComparer : IEqualityComparer<Complex>
        {
            public bool Equals(Complex x, Complex y)
            {
                return (x - y).Magnitude / (x + y).Magnitude < PRECISION;
            }

            public int GetHashCode(Complex obj)
            {
                throw new NotImplementedException();
            }
        }
        protected static IEqualityComparer<Complex> complexComparer = new ComplexComparer();

        public class DftTheoryTestDataProvider : IEnumerable<object[]>
        {
            public IEnumerator<object[]> GetEnumerator()
            {
                return GenerateTestData().GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return GetEnumerator();
            }

            private IEnumerable<object[]> GenerateTestData()
            {
                yield return SingleDataPoint;
                yield return DualDataPoint;
                yield return QuadDataPoint;
                yield return new object[] { FftTestData4096.Signal, FftTestData4096.Spectrum };
            }

            private object[] SingleDataPoint
            {
                get
                {
                    Complex[] data = new Complex[]
                    {
                        new Complex(1, 1)
                    };
                    return new object[] { data, data };
                }
            }

            private object[] DualDataPoint
            {
                get
                {
                    Complex[] signal = new Complex[]
                    {
                        new Complex(1.7, 3.2),
                        new Complex(-2.3, 1.6)
                    };

                    Complex[] spectrum = new Complex[]
                    {
                        new Complex(-0.6, 4.8),
                        new Complex(4, 1.6)
                    };

                    return new object[] { signal, spectrum };
                }
            }

            private object[] QuadDataPoint
            {
                get
                {
                    Complex[] signal = new Complex[]
                    {
                        new Complex(0.3, -2.7),
                        new Complex(-2.2, -1.8),
                        new Complex(3.1, -2.5),
                        new Complex(3.1, 0.7)
                    };

                    Complex[] spectrum = new Complex[]
                    {
                        new Complex(4.3, -6.3),
                        new Complex(-5.3, 5.1),
                        new Complex(2.5, -4.1),
                        new Complex(-0.3, -5.5)
                    };

                    return new object[] { signal, spectrum };
                }
            }
        }

        [Theory]
        [ClassData(typeof(DftTheoryTestDataProvider))]
        public void ArrayDftTheory(Complex[] signalExpected, Complex[] spectrumExpected)
        {
            int windowSize = signalExpected.Length;
            IFftCalculator calculator = CreateCalculator(windowSize);

            Complex[] spectrumActual = new Complex[windowSize];
            calculator.DFT(signalExpected, spectrumActual);
            Assert.Equal(spectrumExpected, spectrumActual, complexComparer);

            Complex[] signalActual = new Complex[windowSize];
            calculator.IDFT(spectrumActual, signalActual);
            Assert.Equal(signalExpected, signalActual, complexComparer);
        }

        [Theory]
        [ClassData(typeof(DftTheoryTestDataProvider))]
        public void SpanDftTheory(Complex[] signalExpected, Complex[] spectrumExpected)
        {
            int windowSize = signalExpected.Length;
            IFftCalculator calculator = CreateCalculator(windowSize);

            Span<double> signalExpectedReal = stackalloc double[windowSize];
            Span<double> signalExpectedImaginary = stackalloc double[windowSize];
            ComplexSpan signalExpectedSpan = new ComplexSpan(signalExpectedReal, signalExpectedImaginary);
            signalExpected.CopyTo(signalExpectedSpan);
            Span<double> spectrumExpectedReal = stackalloc double[windowSize];
            Span<double> spectrumExpectedImaginary = stackalloc double[windowSize];
            ComplexSpan spectrumExpectedSpan = new ComplexSpan(spectrumExpectedReal, spectrumExpectedImaginary);
            spectrumExpected.CopyTo(spectrumExpectedSpan);

            Span<double> signalActualReal = stackalloc double[windowSize];
            Span<double> signalActualImaginary = stackalloc double[windowSize];
            ComplexSpan signalActualSpan = new ComplexSpan(signalActualReal, signalActualImaginary);
            Span<double> spectrumActualReal = stackalloc double[windowSize];
            Span<double> spectrumActualImaginary = stackalloc double[windowSize];
            ComplexSpan spectrumActualSpan = new ComplexSpan(spectrumActualReal, spectrumActualImaginary);

            calculator.DFT(signalExpectedSpan, spectrumActualSpan);
            Assert.Equal(spectrumExpected, spectrumActualSpan.ToArray(), complexComparer);

            calculator.IDFT(spectrumActualSpan, signalActualSpan);
            Assert.Equal(signalExpected, signalActualSpan.ToArray(), complexComparer);
        }
    }
}
