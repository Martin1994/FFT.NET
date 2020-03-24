using System;
using System.Buffers;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;

namespace MartinCl2.Numerics
{
    public sealed class ParallelFftCalculator : IFftCalculator
    {
        private readonly int size;
        private readonly int depth;
        private readonly Complex[] coefficientDFT;
        private readonly Complex[] coefficientIDFT;
        private readonly Complex[] temp;

        public ParallelFftCalculator(int size)
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

        public void IDFT(Complex[] spectrum, Complex[] signal)
        {
            AbstractFFT(spectrum, signal, coefficientIDFT);
            Complex coefficientAdjustment = Complex.One / size;
            for (int i = 0; i < size; i++)
            {
                signal[i] *= coefficientAdjustment;
            }
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

            // We need O(n) extra memory. Make sure the last layer is assigned to the output array.
            Complex[] lowerLayer = depth % 2 == 0 ? to : temp;
            Complex[] upperLayer = depth % 2 == 0 ? temp : to;

            int taskNum = Environment.ProcessorCount;
            Task[] tasks = new Task[taskNum - 1];
            long[] partition = new long[taskNum + 1];
            for (int i = 0; i < taskNum; i++)
            {
                partition[i] = size / taskNum * i;
            }
            partition[taskNum] = size;

            // Assign the initial layer.
            Array.Copy(from, lowerLayer, WindowSize);

            // Calculate the remaining layers.
            // i.e. calculating the ith layer based on the (i + 1)th layer.
            for (int i = 1; i <= depth; i++)
            {
                long p = 0;
                for (int t = 0; t < taskNum - 1; t++)
                {
                    tasks[t] = Task.Run(() =>
                    {
                        long thisP = Interlocked.Increment(ref p);
                        long jFrom = partition[thisP];
                        long jTo = partition[thisP + 1];
                        FFTWorker(upperLayer, lowerLayer, coefficient, i, jFrom, jTo);
                    });
                }
                FFTWorker(upperLayer, lowerLayer, coefficient, i, partition[0], partition[1]);
                Task.WaitAll(tasks);

                Complex[] exchange = upperLayer;
                upperLayer = lowerLayer;
                lowerLayer = exchange;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void FFTWorker(Complex[] upperLayer, Complex[] lowerLayer, Complex[] coefficient, int i, long jFrom, long jTo)
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
