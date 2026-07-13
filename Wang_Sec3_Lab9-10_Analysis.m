% =========================================================================
% Physics 326 Lab 9-10: Fourier Series Analysis
% Author: Leran Wang
% Description:
%   Fourier series synthesis and FFT-based spectrum analysis for square,
%   triangular, sawtooth, and acoustic waveforms. Compares MATLAB-computed
%   FFT coefficients against hardware spectra from SR760 FFT Spectrum
%   Analyzer acquired via oscilloscope CSV exports and GPIB data files.
% =========================================================================

% -------------------------------------------------------------------------
% DATA LOADING — Part D: Function Generator Waveforms (Square & Triangular)
% Each waveform type has 3 repeated measurements (D11/D12/D13, D21/D22/D23)
% CSV columns: col 4 = time, col 5 = voltage (oscilloscope export format)
% -------------------------------------------------------------------------
fin=readtable('D11/F0001CH1.CSV');
tSq=fin{:,4};
vSq=fin{:,5};
fin=readtable('D12/F0002CH1.CSV');
tSq=[tSq,fin{:,4}];           % Concatenate 3 runs as columns
vSq=[vSq,fin{:,5}];
fin=readtable('D13/F0003CH1.CSV');
tSq=[tSq,fin{:,4}];
vSq=[vSq,fin{:,5}];

fin=readtable('D21/F0000CH1.CSV');
tTrig=fin{:,4};
vTrig=fin{:,5};
fin=readtable('D22/F0001CH1.CSV');
tTrig=[tTrig,fin{:,4}];
vTrig=[vTrig,fin{:,5}];
fin=readtable('D23/F0002CH1.CSV');
tTrig=[tTrig,fin{:,4}];
vTrig=[vTrig,fin{:,5}];

% Time window for one full period of the generated waveform [s]
Tf=5.38e-5; Ti=-4.76e-5;

% Preserve originals before any slicing/shifting
tSqOrig=tSq; vSqOrig=vSq;

% Load SR760 spectrum data for square and triangular waves (GPIB export)
% Format: col 1 = frequency [Hz], col 2 = amplitude [V]
fin=load('partD1.txt');
freqD1=fin(:,1); ampSRD1=fin(:,2);
fin=load('partD2.txt');
freqD2=fin(:,1); ampSRD2=fin(:,2);

tTrigOrig=tTrig; vTrigOrig=vTrig;

% -------------------------------------------------------------------------
% SECTION 1 — Fourier Series Synthesis: Square Wave
% Coefficient: A_n = -4/(n*pi) * (-1)^((n+1)/2) for odd n; zero otherwise
% Using complex exponential form: f(t) = sum f_hat_n * exp(i*2*pi*n*v0*t)
% f_hat_n = A_n/2 for n >= 0 (symmetric/even function, so B_n = 0)
% -------------------------------------------------------------------------
fn_Sq=@(n)-4/n/pi*(-1)^((n+1)/2)/2;   % Complex Fourier coefficient for square wave

t=linspace(-.5,.5,500).';              % Time axis: one period T=1s
n=[1;3;5;7;9];                         % Odd harmonics to include
fSq=zeros(size(t,1),size(n,1));        % Preallocate: each column = one n_max
T=1; v0=1/T;                           % Fundamental frequency

% Build partial sums up to each n_max; sum over +/- frequencies for real output
for j=1:size(n,1)
    for i=-n(j):2:n(j)                 % Step by 2: only odd harmonics
        fSq(:,j)=fSq(:,j)+fn_Sq(abs(i))*exp((i*2*pi*v0*t)*1i);
    end
end

figure;
plot(t,fSq);
legend('n:1','3','5','7','9');
xlabel('t/s'); ylabel('amplitude: f(t)'); title('modeled even square wave for 1 period');

% -------------------------------------------------------------------------
% SECTION 2 — Fourier Series Synthesis: Sawtooth Wave
% f(t) = 2t on [-T/2, T/2]; odd function so only B_n terms survive
% B_n = -2/(pi*n), giving f(t) = -2/pi * sum sin(2*pi*n*v0*t)/n
% -------------------------------------------------------------------------
n=[1;3;5];
f2t=zeros(size(t,1),size(n,1));

for j=1:size(n,1)
    for i=1:n(j)
        f2t(:,j)=f2t(:,j)-2/pi/i*sin(2*pi*i*v0*t);
    end
end

figure;
plot(t,f2t);
fst=sawtooth(2*pi*v0*t);               % MATLAB built-in for exact comparison
hold on;
plot(t,fst,'r-.');
xlabel('t/s'); ylabel('amplitude: f(t)'); title('modeled and exact sawtooth wave for 1 period');
legend('n=1','3','5','exact');

% -------------------------------------------------------------------------
% SECTION 3 — Generated Square Wave: Raw Oscilloscope Data (3 runs)
% Time-shift by half period to center waveform for model comparison
% -------------------------------------------------------------------------
figure;
tSq=tSqOrig-Tf/2;                      % Center time axis on waveform midpoint
plot(tSq(:,1),vSq(:,1));
hold on;
plot(tSq(:,2),vSq(:,2),'-.');
plot(tSq(:,3),vSq(:,3),':','LineWidth',2);
xlabel('t/s'); ylabel('amplitude Voltage/V'); title('square waves (generated and theoretical)');

% -------------------------------------------------------------------------
% SECTION 4 — Model vs Generated Square Wave Overlay (n=21)
% Reconstructs Fourier model scaled to match measured voltage amplitude,
% then overlays on the oscilloscope waveform for visual validation
% -------------------------------------------------------------------------
t=tSq(:,1);
fSq=zeros(size(t));
v0=1/(Tf-Ti);                           % Fundamental freq from measured period

for i=-21:2:21                          % Sum harmonics up to n=21
    fSq=fSq+fn_Sq(abs(i))*exp((i*2*pi*v0*t)*1i);
end

% Slice to one period window and scale model amplitude to match data
tSlice=t(t<=Tf & t>=Ti);
fSlice=fSq(t<=Tf & t>=Ti);
scale=mean((max(vSq)-min(vSq))/2);     % Average peak-to-peak / 2 across runs
fSq=fSq*scale;

figure;
plot(tSlice,fSlice);
hold on;

tSq=tSq(:,1);
vSq=vSq(:,1);
vSlice=vSq(tSq<=Tf & tSq>=Ti);
[tSorted,ind]=sort(tSlice);
vSorted=vSlice(ind);
plot(tSlice,vSlice);
xlabel('t/s'); ylabel('Amplitude voltage/V'); title('modeled and generated square wave');
legend('model','generated');

% -------------------------------------------------------------------------
% SECTION 5 — FFT of Square Wave vs SR760 Spectrum
% Computes single-sided FFT from oscilloscope data (run 3),
% normalizes to recover Fourier coefficient amplitudes,
% and overlays against SR760 hardware FFT spectrum
% Note: frequency axis shifted by -46250 Hz to align DC bins between
%       oscilloscope sample rate and SR760 center frequency
% -------------------------------------------------------------------------
tSq=tSqOrig(:,3);
vSq=vSqOrig(:,3);

Fs=1/mean(diff(tSq));                  % Sampling frequency [Hz]
N=length(vSq);
f=(0:N/2-1)*(Fs/N);                    % Single-sided frequency bins
fftSignal=fft(vSq);
amp=abs(fftSignal)/N*2;                % Normalize: factor of 2 for single-sided, /N for amplitude
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
f=f-46250;                             % Frequency alignment correction (DC offset between instruments)
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal');
hold on;
plot(freqD1,ampSRD1);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 6 — Generated Triangular Wave: Raw Oscilloscope Data (3 runs)
% -------------------------------------------------------------------------
figure;
tTrig=tTrigOrig-Tf/2;
plot(tTrig(:,1),vTrig(:,1));
hold on;
plot(tTrig(:,2),vTrig(:,2),'-.');
plot(tTrig(:,3),vTrig(:,3),':','LineWidth',2);
xlabel('t/s'); ylabel('amplitude Voltage/V'); title('Triangle waves (generated and theoretical)');

% -------------------------------------------------------------------------
% SECTION 7 — Model vs Generated Triangular Wave Overlay (n=21)
% Coefficient: A_n = 8/(pi^2 * n^2) for odd n (even symmetric function)
% -------------------------------------------------------------------------
t=tTrig(:,1);
fTrig=zeros(size(t));
v0=1/(Tf-Ti);
fn_Trig=@(n)4/pi^2/n^2;               % Fourier coefficient for triangular wave

for i=-21:2:21
   fTrig=fTrig+fn_Trig(abs(i))*exp((i*2*pi*v0*t)*1i);
end

scale=mean((max(vTrig)-min(vTrig))/2);
fTrig=fTrig*scale;
tSlice=t(t<=Tf & t>=Ti);
fSlice=fTrig(t<=Tf & t>=Ti);

figure;
plot(tSlice,fSlice);
hold on;

vTmp=vTrig(:,1);
vSlice=vTmp(t<=Tf & t>=Ti);
plot(tSlice,vSlice);
xlabel('t/s'); ylabel('Amplitude voltage/V'); title('modeled and generated Triguare wave');
legend('model','generated');

% -------------------------------------------------------------------------
% SECTION 8 — FFT of Triangular Wave vs SR760 Spectrum
% Same pipeline as Section 5; triangular wave has faster harmonic rolloff
% (A_n ~ 1/n^2) so higher harmonics are suppressed more than square wave
% -------------------------------------------------------------------------
tTrig=tTrigOrig(:,3);
vTrig=vTrigOrig(:,3);

Fs=1/mean(diff(tTrig));
N=length(vTrig);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vTrig);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
f=f-46250;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal');
hold on;
plot(freqD2,ampSRD2);
legend('FFT','SR760');

% =========================================================================
% PART E — Acoustic Analysis: Tuning Forks
% Three tuning forks measured: 440 Hz, 293.67 Hz, 493.88 Hz
% Each has 3 repeated oscilloscope recordings + SR760 linear spectrum
% Expected: near-pure sinusoid → single dominant FFT peak at labeled frequency
% =========================================================================

% --- Load tuning fork 440 Hz data ---
fin=readtable('E11/F0002CH1.CSV'); tE1=fin{:,4}; vE1=fin{:,5};
fin=readtable('E12/F0001CH1.CSV'); tE1=[tE1,fin{:,4}]; vE1=[vE1,fin{:,5}];
fin=readtable('E13/F0000CH1.CSV'); tE1=[tE1,fin{:,4}]; vE1=[vE1,fin{:,5}];
tE1Orig=tE1; vE1Orig=vE1;

fin=load('partE1lin.txt');
freqE1=fin(:,1); ampE1=fin(:,2); freqE1Orig=freqE1; ampE1Orig=ampE1;

% -------------------------------------------------------------------------
% SECTION 9 — FFT of 440 Hz Tuning Fork (Linear Scale)
% -------------------------------------------------------------------------
tE1=tE1Orig(:,3);
vE1=vE1Orig(:,3);

Fs=1/mean(diff(tE1));
N=length(vE1);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vE1);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for tune fork of 440 hz');
hold on;
plot(freqE1,ampE1);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 10 — FFT of 440 Hz Tuning Fork (Log/dB Scale)
% Converts MATLAB FFT amplitude to dB using mag2db for comparison with
% SR760 log-scale output; manual frequency axis remapping needed because
% log10(0) is undefined at DC bin, handled by linear interpolation
% -------------------------------------------------------------------------
fin=load('partE1log.txt');
freqE1=fin(:,1); ampE1=fin(:,2);

figure;
flog=log10(f);
flog(1)=0;                              % Handle DC bin (log10(0) = -Inf)
flog(14:end)=flog(14:end)-flog(13);    % Remove offset after bin 13
flog(1:13)=flog(1:13)*flog(14)/flog(13); % Rescale lower bins proportionally
flog=flog*max(freqE1)/max(flog);       % Scale to match SR760 frequency range

ampLog=mag2db(amp);                     % Convert linear amplitude to dB
% Shift dB curve to align median level with SR760 output for overlay
plot(flog,ampLog+median(ampE1)-median(ampLog));
xlabel('Frequency (Hz)','LineWidth',2); ylabel('Amplitude/V');
title('log view FFT of Signal for tune fork of 440 hz');
hold on;
plot(freqE1,ampE1);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 11 — FFT of 293.67 Hz Tuning Fork (Linear Scale)
% -------------------------------------------------------------------------
fin=readtable('E21/F0000CH1.CSV'); tE2=fin{:,4}; vE2=fin{:,5};
fin=readtable('E22/F0001CH1.CSV'); tE2=[tE2,fin{:,4}]; vE2=[vE2,fin{:,5}];
fin=readtable('E23/F0002CH1.CSV'); tE2=[tE2,fin{:,4}]; vE2=[vE2,fin{:,5}];
tE2Orig=tE2; vE2Orig=vE2;

fin=load('partE2lin.txt');
freqE2=fin(:,1); ampE2=fin(:,2);

tE2=tE2Orig(:,3);
vE2=vE2Orig(:,3);

Fs=1/mean(diff(tE2));
N=length(vE2);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vE2);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for tune fork of 293.67 hz');
hold on;
plot(freqE2,ampE2);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 12 — FFT of 493.88 Hz Tuning Fork (Linear Scale)
% -------------------------------------------------------------------------
fin=readtable('E31/F0000CH1.CSV'); tE3=fin{:,4}; vE3=fin{:,5};
fin=readtable('E32/F0001CH1.CSV'); tE3=[tE3,fin{:,4}]; vE3=[vE3,fin{:,5}];
fin=readtable('E33/F0002CH1.CSV'); tE3=[tE3,fin{:,4}]; vE3=[vE3,fin{:,5}];
tE3Orig=tE3; vE3Orig=vE3;

fin=load('partE3lin');
freqE3=fin(:,1); ampE3=fin(:,2);

tE3=tE3Orig(:,3);
vE3=vE3Orig(:,3);

Fs=1/mean(diff(tE3));
N=length(vE3);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vE3);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for tune fork of 493.88 hz');
hold on;
plot(freqE3,ampE3);
legend('FFT','SR760');

% =========================================================================
% PART F — Acoustic Analysis: Human Voice
% Three voice recordings: voice 2 = spoken, voice 1 & 3 = sung single tones
% Expected: sung tones show periodic harmonic peaks; spoken voice shows
%           broadband distribution with less regular harmonic structure
% =========================================================================

% -------------------------------------------------------------------------
% SECTION 13 — FFT of Voice 1 (Sung Tone)
% -------------------------------------------------------------------------
fin=readtable('F11/F0000CH1.CSV'); tF1=fin{:,4}; vF1=fin{:,5};
fin=readtable('F12/F0001CH1.CSV'); tF1=[tF1,fin{:,4}]; vF1=[vF1,fin{:,5}];
fin=readtable('F13/F0002CH1.CSV'); tF1=[tF1,fin{:,4}]; vF1=[vF1,fin{:,5}];
tF1Orig=tF1; vF1Orig=vF1;

fin=load('partF1lin');
freqF1=fin(:,1); ampF1=fin(:,2);

tF1=tF1Orig(:,3);
vF1=vF1Orig(:,3);

Fs=1/mean(diff(tF1));
N=length(vF1);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vF1);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for voice 1');
hold on;
plot(freqF1,ampF1);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 14 — FFT of Voice 2 (Spoken Sound)
% Expect broadband spectrum vs. the harmonic peaks of sung tones
% -------------------------------------------------------------------------
fin=readtable('F21/F0000CH1.CSV'); tF2=fin{:,4}; vF2=fin{:,5};
fin=readtable('F22/F0001CH1.CSV'); tF2=[tF2,fin{:,4}]; vF2=[vF2,fin{:,5}];
fin=readtable('F23/F0003CH1.CSV'); tF2=[tF2,fin{:,4}]; vF2=[vF2,fin{:,5}];
tF2Orig=tF2; vF2Orig=vF2;

fin=load('partF2lin');
freqF2=fin(:,1); ampF2=fin(:,2);

tF2=tF2Orig(:,3);
vF2=vF2Orig(:,3);

Fs=1/mean(diff(tF2));
N=length(vF2);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vF2);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for voice 2');
hold on;
plot(freqF2,ampF2);
legend('FFT','SR760');

% -------------------------------------------------------------------------
% SECTION 15 — FFT of Voice 3 (Sung Tone)
% -------------------------------------------------------------------------
fin=readtable('F31/F0000CH1.CSV'); tF3=fin{:,4}; vF3=fin{:,5};
fin=readtable('F32/F0001CH1.CSV'); tF3=[tF3,fin{:,4}]; vF3=[vF3,fin{:,5}];
fin=readtable('F33/F0002CH1.CSV'); tF3=[tF3,fin{:,4}]; vF3=[vF3,fin{:,5}];
tF3Orig=tF3; vF3Orig=vF3;

fin=load('partF3lin');
freqF3=fin(:,1); ampF3=fin(:,2);

tF3=tF3Orig(:,3);
vF3=vF3Orig(:,3);

Fs=1/mean(diff(tF3));
N=length(vF3);
f=(0:N/2-1)*(Fs/N);
fftSignal=fft(vF3);
amp=abs(fftSignal)/N*2;
f=f(1:floor(N/2));
amp=amp(1:floor(N/2));

figure;
plot(f,amp);
xlabel('Frequency (Hz)'); ylabel('Amplitude/V'); title('FFT of Signal for voice 3');
hold on;
plot(freqF3,ampF3);
legend('FFT','SR760');
