using MartinCl2.Numerics;
using System;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.Intrinsics.X86;

namespace Playground
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Number of cores: {0}", Environment.ProcessorCount);
            Console.WriteLine("Accelerated vectors: {0}", Vector.IsHardwareAccelerated);
            Console.WriteLine("Size of vectors: {0}", Vector<double>.Count);
            Console.WriteLine("Support Avx2: {0}", Avx2.IsSupported);
            Console.WriteLine("Support Avx: {0}", Avx.IsSupported);
            Console.WriteLine("Support Sse: {0}", Sse.IsSupported);
            Console.WriteLine("Support Sse2: {0}", Sse2.IsSupported);
            Console.WriteLine("Support Sse3: {0}", Sse3.IsSupported);
            Console.WriteLine("Support Ssse3: {0}", Ssse3.IsSupported);
            Console.WriteLine("Support Sse41: {0}", Sse41.IsSupported);
            Console.WriteLine("Support Sse42: {0}", Sse42.IsSupported);
            Console.WriteLine("Support Fma: {0}", Fma.IsSupported);

            int windowSize = 4096;

            IFftCalculator[] calculatorsUsingArray = new IFftCalculator[]
            {
                new FftCalculator(windowSize),
                new ParallelFftCalculator(windowSize)
            };

            foreach (var harness in calculatorsUsingArray.Select(calculator => new FftBenchmarkHarness(calculator)))
            {
                TimeSpan elapsedTime = harness.BenchmarkUsingArray();
                Console.WriteLine("{0}: {1}", harness.Name, elapsedTime);
            }

            IFftCalculator[] calculatorsUsingSpan = new IFftCalculator[]
            {
                new SpanFftCalculator(windowSize),
                new VectorTFftCalculator(windowSize),
                new AlignedVectorTFftCalculator(windowSize),
                new IntrinsicsFftCalculator(windowSize),
                new AlignedIntrinsicsFftCalculator(windowSize),
            };

            foreach (var harness in calculatorsUsingSpan.Select(calculator => new FftBenchmarkHarness(calculator)))
            {
                TimeSpan elapsedTime = harness.BenchmarkUsingSpan();
                Console.WriteLine("{0}: {1}", harness.Name, elapsedTime);
            }
        }
    }
}
