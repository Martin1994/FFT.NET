using System;
using System.Numerics;

namespace MartinCl2.Numerics
{
    public ref struct ComplexSpan
    {
        public readonly Span<double> Real;
        public readonly Span<double> Imaginary;

        public ComplexSpan(Span<double> real, Span<double> imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public Complex this[int index]
        {
            get => new Complex(Real[index], Imaginary[index]);
            set
            {
                Real[index] = value.Real;
                Imaginary[index] = value.Imaginary;
            }
        }
    }
}
