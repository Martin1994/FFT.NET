using System;
using System.Numerics;

namespace MartinCl2.Numerics
{
    internal static class ComplexSpanHelper
    {
        public static void CopyTo(this Complex[] array, ComplexSpan span)
        {
            Span<double> realSpan = span.Real;
            Span<double> imaginarySpan = span.Imaginary;
            for (int i = 0; i < realSpan.Length; i++)
            {
                realSpan[i] = array[i].Real;
                imaginarySpan[i] = array[i].Imaginary;
            }
        }

        public static void CopyTo(this ReadOnlyComplexSpan span, Complex[] array)
        {
            ReadOnlySpan<double> realSpan = span.Real;
            ReadOnlySpan<double> imaginarySpan = span.Imaginary;
            for (int i = 0; i < realSpan.Length; i++)
            {
                array[i] = new Complex(realSpan[i], imaginarySpan[i]);
            }
        }

        public static void CopyTo(this ComplexSpan span, Complex[] array)
        {
            Span<double> realSpan = span.Real;
            Span<double> imaginarySpan = span.Imaginary;
            for (int i = 0; i < realSpan.Length; i++)
            {
                array[i] = new Complex(realSpan[i], imaginarySpan[i]);
            }
        }

        public static Complex[] ToArray(this ReadOnlyComplexSpan span)
        {
            Complex[] array = new Complex[span.Real.Length];
            span.CopyTo(array);
            return array;
        }

        public static Complex[] ToArray(this ComplexSpan span)
        {
            Complex[] array = new Complex[span.Real.Length];
            span.CopyTo(array);
            return array;
        }
    }
}
