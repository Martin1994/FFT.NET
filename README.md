# A C# FFT implementation

This project contains 1D FFT implementations with benchmark.

Available implementations are:
* `FftCalculator` - Takes `Complex[]` as parameters. No explicit SIMD optimization.
* `ParallelFftCalculator` - Takes `Complex[]` as parameters. Multi-threaded. No explicit SIMD optimization.
* `SpanFftCalculator` - Takes `Span<double>` as parameters. No explicit SIMD optimization.
* `VectorTFftCalculator` - Takes `Span<double>` as parameters. Uses `Vector<double>` under the hood as SIMD optimization.
* `IntrinsicsFftCalculator` - Takes `Span<double>` as parameters. Uses Intel SIMD intrinsics under the hood as SIMD optimization.

# Example usage

Public interface can be found in `IFftCalculator`.

```csharp
// Some constants
int windowSize = 4096;
double[] audioSource = new double[60 * 44100];
double[] audioOutput = new double[60 * 44100];
VectorFftCalculator fft = new VectorFftCalculator(windowSize);

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

# Benchmark

Benchmark code is available in `./Playground`. The benckmark only tests FFT performance with 4096 window size and real input. Please make sure run the benchmark with Release configuration.

## i7-9700K - 3.60GHz / 8C8T - Windows 10
```
Number of cores: 8
Accelerated vectors: True
Size of vectors: 4
Support Avx2: True
Support Avx: True
Support Sse: True
Support Sse2: True
Support Sse3: True
Support Ssse3: True
Support Sse41: True
Support Sse42: True
Support Fma: True
FftCalculator: 00:00:00.0000225
ParallelFftCalculator: 00:00:00.0000479
SpanFftCalculator: 00:00:00.0000254
VectorTFftCalculator: 00:00:00.0000167
IntrinsicsFftCalculator: 00:00:00.0000169
```

Note: numpy on the same machine took 0.0000419 seconds in average.

## i7-9700K - 3.60GHz / 8C8T - Ubuntu 18.04 on WSL1
```
Number of cores: 8
Accelerated vectors: True
Size of vectors: 4
Support Avx2: True
Support Avx: True
Support Sse: True
Support Sse2: True
Support Sse3: True
Support Ssse3: True
Support Sse41: True
Support Sse42: True
Support Fma: True
FftCalculator: 00:00:00.0000237
ParallelFftCalculator: 00:00:00.0000879
SpanFftCalculator: 00:00:00.0000265
VectorTFftCalculator: 00:00:00.0000175
IntrinsicsFftCalculator: 00:00:00.0000182
```

## i5-4260U - 1.40GHz / 2C4T - Windows 10
```
Number of cores: 8
Accelerated vectors: True
Size of vectors: 4
Support Avx2: True
Support Avx: True
Support Sse: True
Support Sse2: True
Support Sse3: True
Support Ssse3: True
Support Sse41: True
Support Sse42: True
Support Fma: True
FftCalculator: 00:00:00.0000999
ParallelFftCalculator: 00:00:00.0001567
SpanFftCalculator: 00:00:00.0001143
VectorTFftCalculator: 00:00:00.0000646
IntrinsicsFftCalculator: 00:00:00.0000638
```

## BCM2711 (Quad-core Cortex-A72) - 1.65GHz / 4C4T - Fedora 31
```
Number of cores: 4
Accelerated vectors: True
Size of vectors: 2
Support Avx2: False
Support Avx: False
Support Sse: False
Support Sse2: False
Support Sse3: False
Support Ssse3: False
Support Sse41: False
Support Sse42: False
Support Fma: False
FftCalculator: 00:00:00.0001050
ParallelFftCalculator: 00:00:00.0002320
SpanFftCalculator: 00:00:00.0001087
VectorTFftCalculator: 00:00:00.0001237
IntrinsicsFftCalculator: 00:00:00.0001139
```

Note: `IntrinsicsFftCalculator` falled back to an SISD implementation

## Platinum 8175M - 2.50GHz / 1C2T virtualized (EC2 t3.micro) - Amazon Linux
```
Number of cores: 2
Accelerated vectors: True
Size of vectors: 4
Support Avx2: True
Support Avx: True
Support Sse: True
Support Sse2: True
Support Sse3: True
Support Ssse3: True
Support Sse41: True
Support Sse42: True
Support Fma: True
FftCalculator: 00:00:00.0000442
ParallelFftCalculator: 00:00:00.0001234
SpanFftCalculator: 00:00:00.0000519
VectorTFftCalculator: 00:00:00.0000391
IntrinsicsFftCalculator: 00:00:00.0000346
```

## Graviton 2 - 2.50GHz / 1C1T virtualized (EC2 m6g.medium) - Amazon Linux
```
Number of cores: 1
Accelerated vectors: True
Size of vectors: 2
Support Avx2: False
Support Avx: False
Support Sse: False
Support Sse2: False
Support Sse3: False
Support Ssse3: False
Support Sse41: False
Support Sse42: False
Support Fma: False
FftCalculator: 00:00:00.0000491
ParallelFftCalculator: 00:00:00.0000669
SpanFftCalculator: 00:00:00.0000481
VectorTFftCalculator: 00:00:00.0000456
IntrinsicsFftCalculator: 00:00:00.0000490
```

Note: `IntrinsicsFftCalculator` falled back to an SISD implementation
