using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MartinCl2.Numerics
{
    public abstract class AbstractFFTCalculator
    {
        protected readonly long size;
        protected readonly int depth;
        private readonly Complex[] coefficientDFT;
        private readonly Complex[] coefficientIDFT;
        
        public AbstractFFTCalculator(long size)
        {
            if (size <= 0)
            {
                throw new ArgumentException("n must be a positive integer.");
            }
            if ((size & (size - 1)) != 0)
            {
                throw new ArgumentException("n must be a power of 2.");
            }
            this.size = size;
            depth = Log2(size);

            coefficientDFT = new Complex[size];
            Complex angleInterval = -2 * Complex.ImaginaryOne * Math.PI / size;
            for (int i = 0; i < size; i++)
            {
                coefficientDFT[i] = Complex.Exp(i * angleInterval);
            }
            coefficientIDFT = new Complex[size];
            for (int i = 0; i < size; i++)
            {
                coefficientIDFT[i] = 1 / (coefficientDFT[i]);
            }
        }

        public void DFT(Complex[] signal, Complex[] spectrum)
        {
            AbstractFFT(signal, spectrum, coefficientDFT);
        }

        public void IDFT(Complex[] spectrum, Complex[] signal)
        {
            AbstractFFT(spectrum, signal, coefficientIDFT);
            Complex coefficientAdjustment = Complex.One / size;
            for (int i = 0; i < size; i++)
            {
                signal[i] *= coefficientAdjustment;
            }
        }

        protected abstract void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient);
        
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        protected void FFTWorker(Complex[] upperLayer, Complex[] lowerLayer, Complex[] coefficient, int i, long jFrom, long jTo)
        {
            // Copy fields onto the stack
            int depth = this.depth;

            long kPeriod = 1L << i;
            long basePeriod = 1L << (depth - i);
            long basePeriodMask = basePeriod - 1;
            long sizeMask = size - 1;
            for (long j = jFrom; j < jTo; j++)
            {
                long k = j >> (depth - i);
                long baseBit = j & basePeriodMask;
                if (k < kPeriod >> 1)
                {
                    Complex even = lowerLayer[baseBit | (k << depth - i + 1)];
                    Complex odd = lowerLayer[baseBit | (k << depth - i + 1) | basePeriod];
                    Complex oddCoefficient = coefficient[k << depth - i & sizeMask];
                    upperLayer[j] = even + oddCoefficient * odd;
                }
                else
                {
                    long reducedK = k - (kPeriod >> 1);
                    Complex even = lowerLayer[baseBit | (reducedK << depth - i + 1)];
                    Complex odd = lowerLayer[baseBit | (reducedK << depth - i + 1) | basePeriod];
                    Complex oddCoefficient = coefficient[reducedK << depth - i & sizeMask];
                    upperLayer[j] = even - oddCoefficient * odd;
                }
            }
        }

        private int Log2(long n)
        {
            int log = 0;
            while (n > 1)
            {
                n = n >> 1;
                log++;
            }
            return log;
        }
    }
}