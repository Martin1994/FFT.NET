using System;
using System.Diagnostics;
using System.Numerics;

namespace MartinCl2.Numerics
{
    public class FftBenchmarkHarness
    {
        private readonly IFftCalculator calculator;
        private readonly Complex[] testInput = new Complex[WINDOW_SIZE];

        private const int PREWARM_TIMES = 10;
        private const int TEST_TIMES = 100000;
        private const int WINDOW_SIZE = 4096;

        public string Name
        {
            get => calculator.GetType().Name;
        }

        public FftBenchmarkHarness(IFftCalculator calculator)
        {
            this.calculator = calculator;

            Random rnd = new Random();
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                testInput[i] = rnd.NextDouble();
            }
        }

        public TimeSpan BenchmarkUsingArray()
        {
            Complex[] output = new Complex[4096];
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

        public TimeSpan BenchmarkUsingSpan()
        {
            Span<double> outputReal = stackalloc double[4096];
            Span<double> outputImaginary = stackalloc double[4096];
            ComplexSpan output = new ComplexSpan(outputReal, outputImaginary);
            Span<double> inputReal = stackalloc double[4096];
            Span<double> inputImaginary = stackalloc double[4096];
            ComplexSpan input = new ComplexSpan(inputReal, inputImaginary);
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
