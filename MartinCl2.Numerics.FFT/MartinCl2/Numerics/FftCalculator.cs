using System;
using System.Diagnostics;
using System.Numerics;

namespace MartinCl2.Numerics
{
    public class FftCalculator
    {
        private readonly int windowSize;
        private readonly int depth;
        private readonly double[] dftCoefficientReal;
        private readonly double[] dftCoefficientImaginary;
        private readonly double[] idftCoefficientReal;
        private readonly double[] idftCoefficientImaginary;

        public FftCalculator(int windowSize)
        {
            if (windowSize <= 0)
            {
                throw new ArgumentException("n must be a positive integer.");
            }
            if ((windowSize & (windowSize - 1)) != 0)
            {
                throw new ArgumentException("n must be a power of 2.");
            }
            this.windowSize = windowSize;
            depth = Log2(windowSize);

            // Initialize coefficient
            Complex angleInterval = -2 * Complex.ImaginaryOne * Math.PI / windowSize;
            dftCoefficientReal = new double[windowSize];
            dftCoefficientImaginary = new double[windowSize];
            for (int i = 0; i < windowSize; i++)
            {
                Complex coefficient = Complex.Exp(i * angleInterval);
                dftCoefficientReal[i] = coefficient.Real;
                dftCoefficientImaginary[i] = coefficient.Imaginary;
            }

            idftCoefficientReal = new double[windowSize];
            idftCoefficientImaginary = new double[windowSize];
            for (int i = 0; i < windowSize; i++)
            {
                Complex coefficient = 1 / Complex.Exp(i * angleInterval);
                idftCoefficientReal[i] = coefficient.Real;
                idftCoefficientImaginary[i] = coefficient.Imaginary;
            }
        }

        public int WindowSize { get => windowSize; }

        public void DFT(ReadOnlySpan<double> signalReal, ReadOnlySpan<double> signalImaginary, Span<double> spectrumReal, Span<double> spectrumImaginary)
        {
            FFT(new ReadOnlyComplexSpan(signalReal, signalImaginary), new ComplexSpan(spectrumReal, spectrumImaginary), new ComplexSpan(dftCoefficientReal, dftCoefficientImaginary));
        }

        public void IDFT(ReadOnlySpan<double> spectrumReal, ReadOnlySpan<double> spectrumImaginary, Span<double> signalReal, Span<double> signalImaginary)
        {
            FFT(new ReadOnlyComplexSpan(spectrumReal, spectrumImaginary), new ComplexSpan(signalReal, signalImaginary), new ComplexSpan(idftCoefficientReal, idftCoefficientImaginary));
            double coefficientAdjustment = 1D / windowSize;
            for (int i = 0; i < windowSize; i++)
            {
                signalReal[i] *= coefficientAdjustment;
                signalImaginary[i] *= coefficientAdjustment;
            }
        }

        private void FFT(ReadOnlyComplexSpan intput, ComplexSpan output, ComplexSpan dftCoefficient)
        {
            Debug.Assert(intput.Real.Length == windowSize);
            Debug.Assert(intput.Imaginary.Length == windowSize);
            Debug.Assert(output.Real.Length == windowSize);
            Debug.Assert(output.Imaginary.Length == windowSize);

            // Cache variables onto the stack
            int depth = this.depth;

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            // Lower layer will be the final output.
            Span<double> tempReal = stackalloc double[windowSize];
            Span<double> tempImaginary = stackalloc double[windowSize];
            ComplexSpan temp = new ComplexSpan(tempReal, tempImaginary);
            ComplexSpan passInput = depth % 2 == 0 ? output : temp;
            ComplexSpan passOutput = depth % 2 == 0 ? temp : output;

            // Assign the initial pass.
            intput.Real.CopyTo(passInput.Real);
            intput.Imaginary.CopyTo(passInput.Imaginary);

            // Calculate the remaining passes.
            // i.e. calculating the ith pass based on the output of the (i - 1)th path.
            int windowSizeMask = windowSize - 1;
            int vectorCount = Vector<double>.Count;
            Span<double> dftCoefficientReal = dftCoefficient.Real;
            Span<double> dftCoefficientImaginary = dftCoefficient.Imaginary;
            for (int pass = 1; pass <= depth; pass++)
            {
                int halfCapacityOfTopBits = 1 << pass - 1; // 2 ^ (pass - 1)
                int capacityOfBottomBits = 1 << (depth - pass); // 2 ^ (depth - pass)
                int bottomBitsMask = capacityOfBottomBits - 1;
                int topBitsMask = windowSizeMask & ~bottomBitsMask;

                Span<double> passInputReal = passInput.Real;
                Span<double> passInputImaginary = passInput.Imaginary;
                Span<double> passOutputReal = passOutput.Real;
                Span<double> passOutputImaginary = passOutput.Imaginary;
                for (int i = 0; i <= windowSizeMask;) // i will increment in nested loops
                {
                    int topIndex = i >> (depth - pass); // i % 2 ^ (depth - #pass), aka top #pass bits
                    if (capacityOfBottomBits < vectorCount)
                    {
                        if (topIndex < halfCapacityOfTopBits)
                        {
                            int coefficientIndex = topIndex << depth - pass;
                            double oddCoefficientReal = dftCoefficientReal[coefficientIndex];
                            double oddCoefficientImaginary = dftCoefficientImaginary[coefficientIndex];
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                double evenReal = passInputReal[evenIndex];
                                double evenImaginary = passInputImaginary[evenIndex];
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                double oddReal = passInputReal[oddIndex];
                                double oddImaginary = passInputImaginary[oddIndex];
                                passOutputReal[i] = evenReal + oddCoefficientReal * oddReal - oddCoefficientImaginary * oddImaginary;
                                passOutputImaginary[i] = evenImaginary + oddCoefficientReal * oddImaginary + oddCoefficientImaginary * oddReal;
                            }
                        }
                        else
                        {
                            int coefficientIndex = topIndex - halfCapacityOfTopBits << depth - pass;
                            double oddCoefficientReal = dftCoefficientReal[coefficientIndex];
                            double oddCoefficientImaginary = dftCoefficientImaginary[coefficientIndex];
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                double evenReal = passInputReal[evenIndex];
                                double evenImaginary = passInputImaginary[evenIndex];
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                double oddReal = passInputReal[oddIndex];
                                double oddImaginary = passInputImaginary[oddIndex];
                                passOutputReal[i] = evenReal - oddCoefficientReal * oddReal + oddCoefficientImaginary * oddImaginary;
                                passOutputImaginary[i] = evenImaginary - oddCoefficientReal * oddImaginary - oddCoefficientImaginary * oddReal;
                            }
                        }
                    }
                    else
                    {
                        if (topIndex < halfCapacityOfTopBits)
                        {
                            int coefficientIndex = topIndex << depth - pass;
                            Vector<double> oddCoefficientReal = new Vector<double>(dftCoefficientReal[coefficientIndex]);
                            Vector<double> oddCoefficientImaginary = new Vector<double>(dftCoefficientImaginary[coefficientIndex]);
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                Vector<double> evenReal = new Vector<double>(passInputReal.Slice(evenIndex));
                                Vector<double> evenImaginary = new Vector<double>(passInputImaginary.Slice(evenIndex));
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                Vector<double> oddReal = new Vector<double>(passInputReal.Slice(oddIndex));
                                Vector<double> oddImaginary = new Vector<double>(passInputImaginary.Slice(oddIndex));
                                (evenReal + oddCoefficientReal * oddReal - oddCoefficientImaginary * oddImaginary).CopyTo(passOutputReal.Slice(i));
                                (evenImaginary + oddCoefficientReal * oddImaginary + oddCoefficientImaginary * oddReal).CopyTo(passOutputImaginary.Slice(i));
                            }
                        }
                        else
                        {
                            int coefficientIndex = topIndex - halfCapacityOfTopBits << depth - pass;
                            Vector<double> oddCoefficientReal = new Vector<double>(dftCoefficientReal[coefficientIndex]);
                            Vector<double> oddCoefficientImaginary = new Vector<double>(dftCoefficientImaginary[coefficientIndex]);
                            int butterflyIndex = coefficientIndex << 1;
                            for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex += vectorCount, i += vectorCount)
                            {
                                int evenIndex = bottomIndex + butterflyIndex;
                                Vector<double> evenReal = new Vector<double>(passInputReal.Slice(evenIndex));
                                Vector<double> evenImaginary = new Vector<double>(passInputImaginary.Slice(evenIndex));
                                int oddIndex = bottomIndex + butterflyIndex + capacityOfBottomBits;
                                Vector<double> oddReal = new Vector<double>(passInputReal.Slice(oddIndex));
                                Vector<double> oddImaginary = new Vector<double>(passInputImaginary.Slice(oddIndex));
                                (evenReal - oddCoefficientReal * oddReal + oddCoefficientImaginary * oddImaginary).CopyTo(passOutputReal.Slice(i));
                                (evenImaginary - oddCoefficientReal * oddImaginary - oddCoefficientImaginary * oddReal).CopyTo(passOutputImaginary.Slice(i));
                            }
                        }
                    }
                }
                ComplexSpan exchange = passOutput;
                passOutput = passInput;
                passInput = exchange;
            }
        }

        private static int Log2(int n)
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
