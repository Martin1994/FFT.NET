using System.Numerics;

namespace MartinCl2.Numerics
{
    public interface IFftCalculator
    {
        int WindowSize { get; }
        void DFT(Complex[] signal, Complex[] spectrum);
        void IDFT(Complex[] spectrum, Complex[] signal);
        void DFT(ReadOnlyComplexSpan signal, ComplexSpan spectrum);
        void IDFT(ReadOnlyComplexSpan spectrum, ComplexSpan signal);
    }
}
