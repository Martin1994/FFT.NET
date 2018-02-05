using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MartinCl2.Numerics
{
    public sealed class FFTCalculator : AbstractFFTCalculator
    {
        private readonly Complex[] temp;

        public FFTCalculator(long size) : base(size)
        {
            temp = new Complex[size];
        }
        
        protected override void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient)
        {
            Debug.Assert(from != null);
            Debug.Assert(from.LongLength == size);
            Debug.Assert(to != null);
            Debug.Assert(to.LongLength == size);

            // Cache depth onto the stack
            int depth = this.depth;

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            // Lower layer will be the final output.
            Complex[] lowerLayer = depth % 2 == 0 ? to : temp;
            Complex[] upperLayer = depth % 2 == 0 ? temp : to;

            // Assign the initial layer.
            from.CopyTo(lowerLayer, 0);

            // Calculate the remaining layers.
            // i.e. calculating the ith layer based on the (i + 1)th layer.
            for (int i = 1; i <= depth; i++)
            {
                long kPeriod = 1L << i;
                long basePeriod = 1L << (depth - i);
                long basePeriodMask = basePeriod - 1;
                long sizeMask = size - 1;
                for (long j = 0; j <= sizeMask; j++)
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
                Complex[] exchange = upperLayer;
                upperLayer = lowerLayer;
                lowerLayer = exchange;
            }
        }
    }
}
