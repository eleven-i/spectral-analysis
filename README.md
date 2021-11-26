# Spectral Analysis

TypeScript library spectral analysis such as PSD and FFT. For calculating FFT, this package uses [kissfft-js](https://www.npmjs.com/package/kissfft-js). This package extends it to make it easier to use and combines it with windowing functions.

This package contains three public functions `calculateFFT`, `welch` and `spectrogram`.

Each of these functions accept the same inputs, an array of data, the sample rate of the data (frequency), the [window](http://en.wikipedia.org/wiki/Window_function) size, the overlap (50% overlap is default), and the windowing function to use.

The windowing functions available are those found in [fft-windowing-ts](https://www.npmjs.com/package/fft-windowing-ts).

---

## FFT

Raw FFT of data. 50% window overlap is standard.

```TS
const calculateFFT = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
)
```

---

## PSD (Welch's Method)

Implementation of Welch's method of spectral density estimation.

```TS
const welch = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
)
```

---

## Spectrogram

Generates a two dimensional array, and set of corresponding frequencies for plotting a spectogram of PSDs

```TS
const spectrogram = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
)
```
