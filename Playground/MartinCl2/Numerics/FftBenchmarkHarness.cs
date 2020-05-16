using System;
using System.Diagnostics;
using System.Numerics;

namespace MartinCl2.Numerics
{
    public class FftBenchmarkHarness
    {
        private readonly IFftCalculator calculator;
        private readonly Complex[] testInput;
        private int WindowSize { get => calculator.WindowSize; }

        private const int PREWARM_TIMES = 10;
        private const int TEST_TIMES = 10000;
        private const int BIT_ALIGNMENT = 512;

        public string Name
        {
            get => calculator.GetType().Name;
        }

        public FftBenchmarkHarness(IFftCalculator calculator)
        {
            this.calculator = calculator;
            testInput = new Complex[WindowSize];

            Random rnd = new Random();
            for (int i = 0; i < WindowSize; i++)
            {
                testInput[i] = rnd.NextDouble();
            }
        }

        public TimeSpan BenchmarkUsingArray()
        {
            Complex[] output = new Complex[WindowSize];
            for (int i = 0; i < PREWARM_TIMES; i++)
            {
                calculator.DFT(testInput, output);
            }

            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; i < TEST_TIMES; i++)
            {
                calculator.DFT(testInput, output);
            }
            return TimeSpan.FromMilliseconds((double)sw.ElapsedMilliseconds / (double)TEST_TIMES);
        }

        public unsafe TimeSpan BenchmarkUsingSpan()
        {
            double* outputMemory = stackalloc double[WindowSize * 2 + BIT_ALIGNMENT / sizeof(double)];
            ComplexSpan output = new ComplexSpan(outputMemory, WindowSize, BIT_ALIGNMENT);
            double* inputMemory = stackalloc double[WindowSize * 2 + BIT_ALIGNMENT / sizeof(double)];
            ComplexSpan input = new ComplexSpan(inputMemory, WindowSize, BIT_ALIGNMENT);
            testInput.CopyTo(input);

            for (int i = 0; i < PREWARM_TIMES; i++)
            {
                calculator.DFT(input, output);
            }

            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; i < TEST_TIMES; i++)
            {
                calculator.DFT(input, output);
            }
            return TimeSpan.FromMilliseconds((double)sw.ElapsedMilliseconds / (double)TEST_TIMES);
        }
    }
}
