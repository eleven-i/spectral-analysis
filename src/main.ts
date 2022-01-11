import { FFTR } from "kissfft-js";
import {
  hann,
  hamming,
  cosine,
  lanczos,
  gaussian,
  tukey,
  blackman,
  exact_blackman,
  kaiser,
  nuttall,
  blackman_harris,
  blackman_nuttall,
  flat_top,
  WindowFunctionName,
} from "fft-windowing-ts";

const findWindowFunction = (name: WindowFunctionName) => {
  switch (name) {
    case "hann":
      return hann;
    case "hamming":
      return hamming;
    case "cosine":
      return cosine;
    case "lanczos":
      return lanczos;
    case "gaussian":
      return gaussian;
    case "tukey":
      return tukey;
    case "blackman":
      return blackman;
    case "exact_blackman":
      return exact_blackman;
    case "kaiser":
      return kaiser;
    case "nuttall":
      return nuttall;
    case "blackman_harris":
      return blackman_harris;
    case "blackman_nuttall":
      return blackman_nuttall;
    case "flat_top":
      return flat_top;
  }
};

const calculateMagnitude = (complexData: number[]) => {
  const newData = [];
  for (let i = 0; i < complexData.length; i = i + 2) {
    newData.push(
      Math.sqrt(
        complexData[i] * complexData[i] +
          complexData[i + 1] * complexData[i + 1]
      )
    );
  }
  return newData;
};

const fft = (inputData: number[]) => {
  const dataLength = inputData.length;
  const fftr = new FFTR(dataLength);
  const transform = fftr.forward(inputData);
  fftr.dispose();
  const magnitude = calculateMagnitude(transform);
  return magnitude;
};

const calculateFFTFreq = (dataLength: number, sampleRate: number) => {
  const fftfreq = Array.from(
    Array(dataLength),
    (x, i) => i * (sampleRate / dataLength)
  );
  return fftfreq;
};

const calculateWindows = (
  inputData: number[],
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName
) => {
  const wf = findWindowFunction(windowingFunction);
  let overlapFactor = 1 / (1 - overlap);
  [overlap, windowSize, overlapFactor] = roundOverlapAndWindowSize(
    overlap,
    windowSize
  );
  const dataLength = inputData.length;
  if (windowSize > dataLength) {
    throw new Error("Window size must be smaller than data size");
  }
  if (overlapFactor > windowSize / 2) {
    throw new Error("Too much overlap for the window size");
  }
  const numberOfWindows =
    Math.floor((overlapFactor * dataLength) / windowSize) - (overlapFactor - 1);
  const windows = [];
  const stepSize = windowSize / overlapFactor;
  //Run window funciton on raw data
  for (let i = 0; i < numberOfWindows; i++) {
    windows.push(wf(inputData.slice(i * stepSize, i * stepSize + windowSize)));
  }
  return windows;
};

const calculatePSDWindows = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName
) => {
  [overlap, windowSize] = roundOverlapAndWindowSize(windowSize, overlap);
  const windows = calculateWindows(
    inputData,
    windowSize,
    overlap,
    windowingFunction
  );
  const scalingFactor =
    1 /
    (sampleRate *
      hann(Array(windowSize).fill(1))
        .map((element) => element ^ 2)
        .reduce((prev, current) => prev + current));
  //Calculate PSD for each window
  const psdWindows = windows.map((window) =>
    fft(window).map((result) => (result ^ 2) * scalingFactor)
  );
  return psdWindows;
};

/**
 * Raw FFT of data. 50% window overlap is standard.
 *
 * @return  {[number[], number[]]}  [Frequencies, FFT]
 */
const calculateFFT = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
) => {
  [overlap, windowSize] = roundOverlapAndWindowSize(windowSize, overlap);
  const windows = calculateWindows(
    inputData,
    windowSize,
    overlap,
    windowingFunction
  );

  const fftWindows = windows.map((window) => fft(window));
  //Combine windows
  const fftResult = fftWindows.reduce(
    (total, current) => current.map((item, i) => total[i] + item),
    Array(windowSize).fill(0)
  );
  const fftfreq = calculateFFTFreq(windowSize, sampleRate);
  return [fftfreq, fftResult];
};

/**
 * Implementation of Welch's method of spectral density estimation.
 *
 * @return  {[number[], number[]]}  [Frequencies, PSD]
 */
const welch = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
) => {
  [overlap, windowSize] = roundOverlapAndWindowSize(windowSize, overlap);
  const psdWindows = calculatePSDWindows(
    inputData,
    sampleRate,
    windowSize,
    overlap,
    windowingFunction
  );
  //Combine windows
  const psd = psdWindows.reduce(
    (total, current) => current.map((item, i) => total[i] + item),
    Array(windowSize).fill(0)
  );
  const fftfreq = calculateFFTFreq(windowSize, sampleRate);
  return [fftfreq, psd];
};

/**
 * Calcuate window size (rounded)
 *
 * @return  number  window size (rounded)
 */
const roundOverlapAndWindowSize = (
  windowSize: number,
  overlap: number
): [number, number, number] => {
  const overlapFactor = Math.round(1 / (1 - overlap));
  //Rounds down the window size to an even number
  const roundedWindowSize =
    overlapFactor * Math.floor(windowSize / overlapFactor);
  const roundedOverlap = 1 / (1 - overlapFactor);
  return [roundedOverlap, roundedWindowSize, overlapFactor];
};

/**
 * Generates a two dimensional array, and set of corresponding frequencies
 * for plotting a spectogram of PSDs
 *
 * @return  {[number[], number[]]}  [Frequencies, PSDs]
 */
const spectrogram = (
  inputData: number[],
  sampleRate: number,
  windowSize: number,
  overlap = 0.5,
  windowingFunction: WindowFunctionName = "hann"
) => {
  [overlap, windowSize] = roundOverlapAndWindowSize(windowSize, overlap);
  const psdWindows = calculatePSDWindows(
    inputData,
    sampleRate,
    windowSize,
    overlap,
    windowingFunction
  );
  const psdWindowTranspose = psdWindows[0].map((x, i) =>
    psdWindows.map((x) => x[i])
  );
  const fftfreq = calculateFFTFreq(windowSize, sampleRate);
  return [fftfreq, psdWindowTranspose];
};

export { calculateFFT, welch, spectrogram };
