# A C# FFT implementation

This project contains a 1D FFT implementation utilizing `Vector<T>`.

This branch contains the most recommended implementation to be used as a library. To see other FFT implementations with benchmark, please go to the master branch.

# Example usage

```csharp
// Some constants
int windowSize = 4096;
double[] audioSource = new double[60 * 44100];
double[] audioOutput = new double[60 * 44100];
FftCalculator fft = new FftCalculator(windowSize);

// Prepare inputs
ReadOnlySpan<double> signalInputReal = new ReadOnlySpan(audioSource, 123 * windowSize, windowSize);
ReadOnlySpan<double> signalInputImaginary = stackalloc double[windowSize];
ReadOnlyComplexSpan signalInput = new ReadOnlyComplexSpan(signalInputReal, signalInputImaginary);

// Prepare outputs
Span<double> spectrumReal = stackalloc double[windowSize];
Span<double> spectrumImaginary = stackalloc double[windowSize];
ComplexSpan specturm = new ComplexSpan(spectrumReal, spectrumImaginary);

// Perform FFT. SignalInput is the input and spectrum is the output.
fft.DFT(signalInput, spectrum);

// Do some magic on the spectrum
// ...

// Prepare outputs
Span<double> signalOutputReal = new Span(audioOutput, 123 * windowSize, windowSize);
Span<double> signalOutputImaginary = stackalloc double[windowSize];
ComplexSpan signalOutput = new ComplexSpan(signalOutputReal, signalOutputImaginary);

// Perform IFFT. Spectrum is the input and signalOutput is the output.
fft.IDFT(spectrum, signalOutput);
```
