using System;
using System.Numerics;

namespace MartinCl2.Numerics
{
    public ref struct ReadOnlyComplexSpan
    {
        public readonly ReadOnlySpan<double> Real;
        public readonly ReadOnlySpan<double> Imaginary;

        public ReadOnlyComplexSpan(ReadOnlySpan<double> real, ReadOnlySpan<double> imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public Complex this[int index]
        {
            get => new Complex(Real[index], Imaginary[index]);
        }

        public static implicit operator ReadOnlyComplexSpan(ComplexSpan rhs)
        {
            return new ReadOnlyComplexSpan(rhs.Real, rhs.Imaginary);
        }
    }
}