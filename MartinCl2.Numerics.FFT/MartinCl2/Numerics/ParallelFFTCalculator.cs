using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;

namespace MartinCl2.Numerics
{
    public sealed class ParallelFFTCalculator : AbstractFFTCalculator
    {
        protected readonly Complex[] temp;
        public ParallelFFTCalculator(long size) : base(size)
        {
            temp = new Complex[size];
        }

        protected override void AbstractFFT(Complex[] from, Complex[] to, Complex[] coefficient)
        {
            Debug.Assert(from != null);
            Debug.Assert(from.Length == size);
            Debug.Assert(to != null);
            Debug.Assert(to.Length == size);

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
            from.CopyTo(lowerLayer, 0);

            // Calculate the remaining layers.
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
    }
}
