# Signal Spectrum Analysis — Fourier Series & FFT Validation

## Overview
MATLAB-based signal analysis pipeline comparing theoretical Fourier 
series models against hardware FFT measurements from an SR760 FFT 
Spectrum Analyzer (acquired via GPIB/LabVIEW), validated across 
multiple waveform types and acoustic signals.

## Equipment
- Stanford Research Systems SR760 FFT Spectrum Analyzer (GPIB)
- Tektronix Oscilloscope
- Function Generator (square/triangular/sawtooth waves)
- Tuning forks: 440 Hz, 293.67 Hz, 493.88 Hz
- Microphone (acoustic signal acquisition)

## What's in this repo
- fourier_synthesis.m  — Fourier series reconstruction up to n=21 
  for square, triangular, sawtooth waves
- fft_analysis.m       — FFT computation from oscilloscope waveform 
  data, normalized to Fourier coefficients
- sr760_compare.m      — Overlay plots of MATLAB FFT vs SR760 hardware 
  spectra (linear + dB scale)
- tuning_fork_analysis.m — Acoustic frequency analysis for 3 tuning forks
- voice_analysis.m     — Harmonic content comparison: sung tone vs spoken voice

## Key Results
- MATLAB FFT and SR760 spectra consistent within measurement uncertainty 
  for square and triangular waves up to ~250 kHz
- Identified spurious even-harmonic peaks in log-scale SR760 output; 
  attributed to non-ideal waveform generation rather than Fourier model error
- Tuning fork measurements confirmed fundamental frequencies within 
  ~0.4% of labeled values (e.g., 441.4 Hz vs 440 Hz nominal)
- Demonstrated harmonic richness difference between sung tones 
  (periodic harmonic peaks) and spoken voice (broadband distribution)
