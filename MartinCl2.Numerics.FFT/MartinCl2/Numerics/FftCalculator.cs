using System;
using System.Buffers;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MartinCl2.Numerics
{
    public class FftCalculator : IFftCalculator
    {
       private readonly int size;
        private readonly int depth;
        private readonly Complex[] coefficientDFT;
        private readonly Complex[] coefficientIDFT;
        private readonly Complex[] temp;

        public FftCalculator(int size)
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

            temp = new Complex[size];
        }

        public int WindowSize { get => size; }

        public void DFT(Complex[] signal, Complex[] spectrum)
        {
            AbstractFFT(signal, spectrum, coefficientDFT);
        }

        public void DFT(ReadOnlyComplexSpan signal, ComplexSpan spectrum)
        {
            Complex[] signalArray = ArrayPool<Complex>.Shared.Rent(size);
            signal.CopyTo(signalArray);
            Complex[] spectrumArray = ArrayPool<Complex>.Shared.Rent(size);

            DFT(signalArray, spectrumArray);

            spectrumArray.CopyTo(spectrum);
            ArrayPool<Complex>.Shared.Return(signalArray);
            ArrayPool<Complex>.Shared.Return(spectrumArray);
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

        public void IDFT(ReadOnlyComplexSpan spectrum, ComplexSpan signal)
        {
            Complex[] signalArray = ArrayPool<Complex>.Shared.Rent(size);
            Complex[] spectrumArray = ArrayPool<Complex>.Shared.Rent(size);
            spectrum.CopyTo(spectrumArray);

            IDFT(spectrumArray, signalArray);

            signalArray.CopyTo(signal);
            ArrayPool<Complex>.Shared.Return(signalArray);
            ArrayPool<Complex>.Shared.Return(spectrumArray);
        }

        private void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient)
        {
            Debug.Assert(from != null);
            Debug.Assert(from.Length >= size);
            Debug.Assert(to != null);
            Debug.Assert(to.Length >= size);

            // Cache depth onto the stack
            int depth = this.depth;

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            // Lower layer will be the final output.
            Complex[] passInput = depth % 2 == 0 ? to : temp;
            Complex[] passOutput = depth % 2 == 0 ? temp : to;

            // Assign the initial layer.
            Array.Copy(from, passInput, WindowSize);

            // Calculate the remaining passes.
            // i.e. calculating the ith pass based on the output of the (i - 1)th path.
            int windowSizeMask = size - 1;
            for (int pass = 1; pass <= depth; pass++)
            {
                int halfCapacityOfTopBits = 1 << pass - 1; // 2 ^ (pass - 1)
                int halfCapacityOfTopBitsMask = halfCapacityOfTopBits - 1;
                int capacityOfBottomBits = 1 << (depth - pass); // 2 ^ (depth - pass)
                int bottomBitsMask = capacityOfBottomBits - 1;
                int topBitsMask = windowSizeMask & ~bottomBitsMask;
                for (int i = 0; i <= windowSizeMask;) // i will increment in nested loops
                {
                    int topIndex = i >> (depth - pass); // i % 2 ^ (depth - #pass), aka top #pass bits
                    if (topIndex < halfCapacityOfTopBits)
                    {
                        Complex oddCoefficient = coefficient[topIndex << depth - pass];
                        int butterflyIndex = topIndex << depth - pass + 1;
                        for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                        {
                            Complex even = passInput[bottomIndex + butterflyIndex];
                            Complex odd = passInput[bottomIndex + butterflyIndex + capacityOfBottomBits];
                            passOutput[i] = even + oddCoefficient * odd;
                        }
                    }
                    else
                    {
                        Complex oddCoefficient = coefficient[topIndex - halfCapacityOfTopBits << depth - pass];
                        int butterflyIndex = topIndex - halfCapacityOfTopBits << depth - pass + 1;
                        for (int bottomIndex = 0; bottomIndex < capacityOfBottomBits; bottomIndex++, i++)
                        {
                            Complex even = passInput[bottomIndex + butterflyIndex];
                            Complex odd = passInput[bottomIndex + butterflyIndex + capacityOfBottomBits];
                            passOutput[i] = even - oddCoefficient * odd;
                        }
                    }
                }
                Complex[] exchange = passOutput;
                passOutput = passInput;
                passInput = exchange;
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
