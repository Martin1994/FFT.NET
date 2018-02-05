# A C# FFT implementation

This project contains a simple FFT implementation by C#.

# Example

```csharp
long size = 4096;

Complex[] signal = new Complex[size];
Complex[] spectrum = new Complex[size];

FFTCalculator fft = new FFTCalculator(size);

// Do an FFT: spectrum = FFT(signal)
fft.DFT(signal, spectrum);

// Do an IFFT: signal = IFFT(spectrum)
fft.IDFT(spectrum, signal);
```
