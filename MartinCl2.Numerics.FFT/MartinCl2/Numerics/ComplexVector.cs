using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace MartinCl2.Numerics
{
    public struct ComplexVector
    {
        public readonly Vector<double> Real;
        public readonly Vector<double> Imaginary;

        public static int Count { get { return Vector<double>.Count; } }

        public Complex this[int i] { get { return new Complex(Real[i], Imaginary[i]); } }

        public ComplexVector(Vector<double> real, Vector<double> imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public ComplexVector(double[] real, double[] imaginary)
        {
            Real = new Vector<double>(real);
            Imaginary = new Vector<double>(imaginary);
        }

        public ComplexVector(double[] real, double[] imaginary, int index)
        {
            Real = new Vector<double>(real, index);
            Imaginary = new Vector<double>(imaginary, index);
        }

        public static ComplexVector operator +(ComplexVector left, ComplexVector right)
        {
            Vector<double> real = left.Real + right.Real;
            Vector<double> imaginary = left.Imaginary + right.Imaginary;
            return new ComplexVector(real, imaginary);
        }

        public static ComplexVector operator -(ComplexVector left, ComplexVector right)
        {
            Vector<double> real = left.Real - right.Real;
            Vector<double> imaginary = left.Imaginary - right.Imaginary;
            return new ComplexVector(real, imaginary);
        }

        public static ComplexVector operator *(ComplexVector left, ComplexVector right)
        {
            Vector<double> real = left.Real * right.Real - left.Imaginary * right.Imaginary;
            Vector<double> imaginary = left.Real * right.Imaginary + left.Imaginary * right.Real;
            return new ComplexVector(real, imaginary);
        }

        public static ComplexVector operator /(ComplexVector left, ComplexVector right)
        {
            Vector<double> denominator = left.Imaginary * left.Imaginary + right.Imaginary * right.Imaginary;
            Vector<double> real = (left.Real * right.Real + left.Imaginary * right.Imaginary) / denominator;
            Vector<double> imaginary = (left.Imaginary * right.Real - left.Real * right.Imaginary) / denominator;
            return new ComplexVector(real, imaginary);
        }
    }
}
