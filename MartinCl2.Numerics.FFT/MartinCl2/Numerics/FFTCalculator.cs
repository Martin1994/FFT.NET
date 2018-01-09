using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MartinCl2.Numerics
{
    public class FFTCalculator
    {
        private readonly long size;
        private readonly int depth;
        private readonly Complex[] coefficientDFT;
        private readonly Complex[] coefficientIDFT;
        private readonly Complex[] temp;
        
        public FFTCalculator(long size)
        {
            if (size <= 0) {
                throw new ArgumentException("n must be a positive integer.");
            }
            if ((size & (size - 1)) != 0) {
                throw new ArgumentException("n must be a power of 2.");
            }
            this.size = size;
            depth = Log2(size);
            temp = new Complex[size];
            coefficientDFT = new Complex[size];
            Complex angleInterval = -2 * Complex.ImaginaryOne * Math.PI / size;
            for (int i = 0; i < size; i++) {
                coefficientDFT[i] = Complex.Exp(i * angleInterval);
            }
            coefficientIDFT = new Complex[size];
            for (int i = 0; i < size; i++) {
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
            for (int i = 0; i < size; i++) {
                signal[i] *= coefficientAdjustment;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient) {
            Debug.Assert(from != null);
            Debug.Assert(from.LongLength == size);
            Debug.Assert(to != null);
            Debug.Assert(to.LongLength == size);

            // Cache depth on stack
            int depth = this.depth;

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            Complex[] lowerLayer = depth % 2 == 0 ? to : temp;
            Complex[] upperLayer = depth % 2 == 0 ? temp : to;

            // Assign the initial layer.
            from.CopyTo(lowerLayer, 0);

            // Calculate the remaining layers.
            for (int i = 1; i <= depth; i++) {
                int kPeriod = 1 << i;
                int basePeriod = 1 << (depth - i);
                int basePeriodMask = basePeriod - 1;
                long sizeMask = size - 1;
                for (int j = 0; j <= sizeMask; j++) {
                    int k = j >> (depth - i);
                    int baseBit = j & basePeriodMask;
                    if (k < kPeriod >> 1) {
                        Complex even = lowerLayer[baseBit | (k << depth - i + 1)];
                        Complex odd = lowerLayer[baseBit | (k << depth - i + 1) | basePeriod];
                        Complex oddCoefficient = coefficient[(long)k << depth - i & sizeMask];
                        upperLayer[j] = even + oddCoefficient * odd;
                    } else {
                        int reducedK = k - (kPeriod >> 1);
                        Complex even = lowerLayer[baseBit | (reducedK << depth - i + 1)];
                        Complex odd = lowerLayer[baseBit | (reducedK << depth - i + 1) | basePeriod];
                        Complex oddCoefficient = coefficient[(long)reducedK << depth - i & sizeMask];
                        upperLayer[j] = even - oddCoefficient * odd;
                    }
                }
                Complex[] exchange = upperLayer;
                upperLayer = lowerLayer;
                lowerLayer = exchange;
            }
        }

        private int Log2(long n) {
            int log = 0;
            while (n > 1) {
                n = n >> 1;
                log++;
            }
            return log;
        }
    }
}
