using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MartinCl2.Numerics
{
    public sealed class VectorizedFFTCalculator : AbstractFFTCalculator
    {
        private readonly ComplexVector[] tempA;
        private readonly ComplexVector[] tempB;

        public VectorizedFFTCalculator(long size) : base(size)
        {
            int vectorSize = Vector<double>.Count;
            if (size % vectorSize != 0)
            {
                throw new PlatformNotSupportedException();
            }

            tempA = new ComplexVector[size / vectorSize];
            tempB = new ComplexVector[size / vectorSize];
        }
        
        protected override void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient)
        {
            Debug.Assert(from != null);
            Debug.Assert(from.LongLength == size);
            Debug.Assert(to != null);
            Debug.Assert(to.LongLength == size);

            int vectorSize = Vector<double>.Count;

            // Cache depth on stack
            int depth = this.depth;

            // Assign buffers
            ComplexVector[] upperLayer = tempA;
            ComplexVector[] lowerLayer = tempB;
            double[] evenReal = new double[vectorSize];
            double[] evenImaginary = new double[vectorSize];
            double[] oddReal = new double[vectorSize];
            double[] oddImaginary = new double[vectorSize];
            double[] coefficientReal = new double[vectorSize];
            double[] coefficientImaginary = new double[vectorSize];

            // Assign the initial layer.
            for (int i = 0, v = 0; i < size; i += vectorSize, v++)
            {
                for (int j = 0; j < vectorSize; j++)
                {
                    evenReal[j] = from[i + j].Real;
                    evenImaginary[j] = from[i + j].Imaginary;
                }
                lowerLayer[v] = new ComplexVector(evenReal, evenImaginary);
            }

            // Calculate the remaining layers.
            for (int i = 1; i <= depth; i++)
            {
                long kPeriod = 1L << i;
                long basePeriod = 1L << (depth - i);
                long basePeriodMask = basePeriod - 1;
                long sizeMask = size - 1;
                int v = 0;
                int arrayIndex = 0;
                for (long j = 0; j <= sizeMask; j++)
                {
                    long k = j >> (depth - i);
                    long baseBit = j & basePeriodMask;
                    if (k < kPeriod >> 1)
                    {
                        // Even
                        long evenIndex = baseBit | (k << depth - i + 1);
                        CopyComplexToBuffer(lowerLayer, evenReal, evenImaginary, vectorSize, evenIndex, v);

                        // Odd
                        long oddIndex = baseBit | (k << depth - i + 1) | basePeriod;
                        CopyComplexToBuffer(lowerLayer, oddReal, oddImaginary, vectorSize, oddIndex, v);

                        // Coefficient
                        long coefficientIndex = k << depth - i & sizeMask;
                        coefficientReal[v] = coefficient[coefficientIndex].Real;
                        coefficientImaginary[v] = coefficient[coefficientIndex].Imaginary;
                    }
                    else
                    {
                        long reducedK = k - (kPeriod >> 1);
                        // Even
                        long evenIndex = baseBit | (reducedK << depth - i + 1);
                        CopyComplexToBuffer(lowerLayer, evenReal, evenImaginary, vectorSize, evenIndex, v);

                        // Odd
                        long oddIndex = baseBit | (reducedK << depth - i + 1) | basePeriod;
                        CopyComplexToBuffer(lowerLayer, oddReal, oddImaginary, vectorSize, oddIndex, v);

                        // Coefficient
                        long coefficientIndex = reducedK << depth - i & sizeMask;
                        coefficientReal[v] = -coefficient[coefficientIndex].Real;
                        coefficientImaginary[v] = -coefficient[coefficientIndex].Imaginary;
                    }

                    // Write buffer to vectors when needed
                    v++;
                    if (v >= vectorSize)
                    {
                        v = 0;
                        ComplexVector even = new ComplexVector(evenReal, evenImaginary);
                        ComplexVector odd = new ComplexVector(oddReal, oddImaginary);
                        ComplexVector oddCoefficient = new ComplexVector(coefficientReal, coefficientImaginary);
                        upperLayer[arrayIndex] = even + oddCoefficient * odd;
                        arrayIndex++;
                    }
                }
                ComplexVector[] exchange = upperLayer;
                upperLayer = lowerLayer;
                lowerLayer = exchange;
            }

            // Copy result back to array
            for (int i = 0, v = vectorSize, arrayIndex = -1; i < size; i++, v++)
            {
                if (v == vectorSize)
                {
                    arrayIndex++;
                    v = 0;
                }
                to[i] = lowerLayer[arrayIndex][v];
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void CopyComplexToBuffer(ComplexVector[] from, double[] realTo, double[] imaginaryTo, long vectorSize, long fromIndex, int v)
        {
            long fromVectorIndex = fromIndex % vectorSize;
            long fromArrayIndex = fromIndex / vectorSize;

            Complex fromValue = from[fromArrayIndex][(int)fromVectorIndex];
            realTo[v] = fromValue.Real;
            imaginaryTo[v] = fromValue.Imaginary;
        }
    }
}
