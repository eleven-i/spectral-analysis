import { WindowFunctionName } from "fft-windowing-ts";
/**
 * Raw FFT of data. 50% window overlap is standard.
 *
 * @return  {[number[], number[]]}  [Frequencies, FFT]
 */
declare const calculateFFT: (inputData: number[], sampleRate: number, windowSize: number, overlap?: number, windowingFunction?: WindowFunctionName) => any[][];
/**
 * Implementation of Welch's method of spectral density estimation.
 *
 * @return  {[number[], number[]]}  [Frequencies, PSD]
 */
declare const welch: (inputData: number[], sampleRate: number, windowSize: number, overlap?: number, windowingFunction?: WindowFunctionName) => any[][];
/**
 * Generates a two dimensional array, and set of corresponding frequencies
 * for plotting a spectogram of PSDs
 *
 * @return  {[number[], number[]]}  [Frequencies, PSDs]
 */
declare const spectrogram: (inputData: number[], sampleRate: number, windowSize: number, overlap?: number, windowingFunction?: WindowFunctionName) => (number[] | number[][])[];
export { calculateFFT, welch, spectrogram };
