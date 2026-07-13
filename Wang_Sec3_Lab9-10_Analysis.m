fin=readtable('D11/F0001CH1.CSV');
tSq=fin{:,4};
vSq=fin{:,5};
fin=readtable('D12/F0002CH1.CSV');
tSq=[tSq,fin{:,4}];
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
Tf=5.38e-5;Ti=-4.76e-5;
tSqOrig=tSq;vSqOrig=vSq;
fin=load('partD1.txt');
freqD1=fin(:,1);ampSRD1=fin(:,2);
fin=load('partD2.txt');
freqD2=fin(:,1);ampSRD2=fin(:,2);
tTrigOrig=tTrig;vTrigOrig=vTrig;
%%
fn_Sq=@(n)-4/n/pi*(-1)^((n+1)/2)/2;
t=linspace(-.5,.5,500).';
n=[1;3;5;7;9];
fSq=zeros(size(t,1),size(n,1));
T=1;
v0=1/T;
for j=1:size(n,1)
    for i=-n(j):2:n(j)
        fSq(:,j)=fSq(:,j)+fn_Sq(abs(i))*exp((i*2*pi*v0*t)*1i);
    end
end
figure;
plot(t,fSq);
legend('n:1','3','5','7','9');
xlabel('t/s');ylabel('amplitude: f(t)');title('modeled even square wave for 1 period');
%%
n=[1;3;5];
f2t=zeros(size(t,1),size(n,1));
for j=1:size(n,1)
    for i=1:n(j)
        f2t(:,j)=f2t(:,j)-2/pi/i*sin(2*pi*i*v0*t);
    end
end
figure;
plot(t,f2t);
fst=sawtooth(2*pi*v0*t);
hold on;
plot(t,fst,'r-.');
xlabel('t/s');ylabel('amplitude: f(t)');title('modeled and exact sawtooth wave for 1 period');
legend('n=1','3','5','exact');
%%
figure;
tSq=tSqOrig-Tf/2;
plot(tSq(:,1),vSq(:,1));
hold on;
plot(tSq(:,2),vSq(:,2),'-.');
plot(tSq(:,3),vSq(:,3),':','LineWidth',2);
xlabel('t/s');ylabel('amplitude Voltage/V');title('square waves (generated and theoretical)');
%%

t=tSq(:,1);
fSq=zeros(size(t));
v0=1/(Tf-Ti);
for i=-21:2:21
    fSq=fSq+fn_Sq(abs(i))*exp((i*2*pi*v0*t)*1i);
end

tSlice=t(t<=Tf & t>=Ti);
fSlice=fSq(t<=Tf & t>=Ti);
scale=mean((max(vSq)-min(vSq))/2);
fSq=fSq*scale;
figure;
plot(tSlice,fSlice);
hold on;
%%
tSq=tSq(:,1);
vSq=vSq(:,1);
vSlice=vSq(tSq<=Tf & tSq>=Ti);
[tSorted,ind]=sort(tSlice);
vSorted=vSlice(ind);
%figure;
plot(tSlice,vSlice);
xlabel('t/s');ylabel('Amplitude voltage/V');title('modeled and generated square wave');
legend('model','generated');
%%

tSq=tSqOrig(:,3);

vSq=vSqOrig(:,3);
Fs=1/mean(diff(tSq));
N=length(vSq);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vSq);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
f=f-46250;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal');
hold on;
plot(freqD1,ampSRD1);
legend('FFT','SR760');
%%
figure;
tTrig=tTrigOrig-Tf/2;
plot(tTrig(:,1),vTrig(:,1));
hold on;
plot(tTrig(:,2),vTrig(:,2),'-.');
plot(tTrig(:,3),vTrig(:,3),':','LineWidth',2);
xlabel('t/s');ylabel('amplitude Voltage/V');title('Triangle waves (generated and theoretical)');
%%
t=tTrig(:,1);
fTrig=zeros(size(t));
v0=1/(Tf-Ti);
fn_Trig=@(n)4/pi^2/n^2;
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
%%
vTmp=vTrig(:,1);
vSlice=vTmp(t<=Tf & t>=Ti);
%figure;
plot(tSlice,vSlice);
xlabel('t/s');ylabel('Amplitude voltage/V');title('modeled and generated Triguare wave');
legend('model','generated');
%%
tTrig=tTrigOrig(:,3);
vTrig=vTrigOrig(:,3);
Fs=1/mean(diff(tTrig));
N=length(vTrig);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vTrig);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
f=f-46250;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal');
hold on;
plot(freqD2,ampSRD2);
legend('FFT','SR760');
%%
fin=readtable('E11/F0002CH1.CSV');
tE1=fin{:,4};
vE1=fin{:,5};
fin=readtable('E12/F0001CH1.CSV');
tE1=[tE1,fin{:,4}];
vE1=[vE1,fin{:,5}];
fin=readtable('E13/F0000CH1.CSV');
tE1=[tE1,fin{:,4}];
vE1=[vE1,fin{:,5}];
tE1Orig=tE1;vE1Orig=vE1;
fin=load('partE1lin.txt');
freqE1=fin(:,1);ampE1=fin(:,2);freqE1Orig=freqE1;ampE1Orig=ampE1;
%%
tE1=tE1Orig(:,3);
vE1=vE1Orig(:,3);
Fs=1/mean(diff(tE1));
N=length(vE1);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vE1);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for tune fork of 440 hz');
hold on;
plot(freqE1,ampE1);
legend('FFT','SR760');
%%
fin=load('partE1log.txt');
freqE1=fin(:,1);ampE1=fin(:,2);
figure;
flog=log10(f);
flog(1)=0;flog(14:end)=flog(14:end)-flog(13);
flog(1:13)=flog(1:13)*flog(14)/flog(13);

flog=flog*max(freqE1)/max(flog);
ampLog=mag2db(amp);
%ampLog=median(ampE1)/median(ampLog)*ampLog;
%ampLog=(max(ampE1)-min(ampE1))/(max(ampLog)-min(ampLog))*ampLog;
plot(flog,ampLog+median(ampE1)-median(ampLog));
xlabel('Frequency (Hz)','LineWidth',2);
ylabel('Amplitude/V');
title('log view FFT of Signal for tune fork of 440 hz');
hold on;
plot(freqE1,ampE1);
legend('FFT','SR760');
%%
fin=readtable('E21/F0000CH1.CSV');
tE2=fin{:,4};
vE2=fin{:,5};
fin=readtable('E22/F0001CH1.CSV');
tE2=[tE2,fin{:,4}];
vE2=[vE2,fin{:,5}];
fin=readtable('E23/F0002CH1.CSV');
tE2=[tE2,fin{:,4}];
vE2=[vE2,fin{:,5}];
tE2Orig=tE2;vE2Orig=vE2;
fin=load('partE2lin.txt');
freqE2=fin(:,1);ampE2=fin(:,2);freqE2Orig=freqE2;ampE2Orig=ampE2;
%%
tE2=tE2Orig(:,3);
vE2=vE2Orig(:,3);
Fs=1/mean(diff(tE2));
N=length(vE2);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vE2);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for tune fork of 293.67 hz');
hold on;
plot(freqE2,ampE2);
legend('FFT','SR760');

%%
fin=readtable('E31/F0000CH1.CSV');
tE3=fin{:,4};
vE3=fin{:,5};
fin=readtable('E32/F0001CH1.CSV');
tE3=[tE3,fin{:,4}];
vE3=[vE3,fin{:,5}];
fin=readtable('E33/F0002CH1.CSV');
tE3=[tE3,fin{:,4}];
vE3=[vE3,fin{:,5}];
tE3Orig=tE3;vE3Orig=vE3;
fin=load('partE3lin');
freqE3=fin(:,1);ampE3=fin(:,2);
%%
tE3=tE3Orig(:,3);
vE3=vE3Orig(:,3);
Fs=1/mean(diff(tE3));
N=length(vE3);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vE3);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for tune fork of 493.88 hz');
hold on;
plot(freqE3,ampE3);
legend('FFT','SR760');

%%
fin=readtable('F11/F0000CH1.CSV');
tF1=fin{:,4};
vF1=fin{:,5};
fin=readtable('F12/F0001CH1.CSV');
tF1=[tF1,fin{:,4}];
vF1=[vF1,fin{:,5}];
fin=readtable('F13/F0002CH1.CSV');
tF1=[tF1,fin{:,4}];
vF1=[vF1,fin{:,5}];
tF1Orig=tF1;vF1Orig=vF1;
fin=load('partF1lin');
freqF1=fin(:,1);ampF1=fin(:,2);
%%
tF1=tF1Orig(:,3);
vF1=vF1Orig(:,3);
Fs=1/mean(diff(tF1));
N=length(vF1);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vF1);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for voice 1');
hold on;
plot(freqF1,ampF1);
legend('FFT','SR760');

%%
fin=readtable('F21/F0000CH1.CSV');
tF2=fin{:,4};
vF2=fin{:,5};
fin=readtable('F22/F0001CH1.CSV');
tF2=[tF2,fin{:,4}];
vF2=[vF2,fin{:,5}];
fin=readtable('F23/F0003CH1.CSV');
tF2=[tF2,fin{:,4}];
vF2=[vF2,fin{:,5}];
tF2Orig=tF2;vF2Orig=vF2;
fin=load('partF2lin');
freqF2=fin(:,1);ampF2=fin(:,2);
%%
tF2=tF2Orig(:,3);
vF2=vF2Orig(:,3);
Fs=1/mean(diff(tF2));
N=length(vF2);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vF2);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for voice 2');
hold on;
plot(freqF2,ampF2);
legend('FFT','SR760');

%%
fin=readtable('F31/F0000CH1.CSV');
tF3=fin{:,4};
vF3=fin{:,5};
fin=readtable('F32/F0001CH1.CSV');
tF3=[tF3,fin{:,4}];
vF3=[vF3,fin{:,5}];
fin=readtable('F33/F0002CH1.CSV');
tF3=[tF3,fin{:,4}];
vF3=[vF3,fin{:,5}];
tF3Orig=tF3;vF3Orig=vF3;
fin=load('partF3lin');
freqF3=fin(:,1);ampF3=fin(:,2);
%%
tF3=tF3Orig(:,3);
vF3=vF3Orig(:,3);
Fs=1/mean(diff(tF3));
N=length(vF3);
f=(0:N/2-1)*(Fs/N); % Frequency bins
fftSignal=fft(vF3);
amp=abs(fftSignal)/N*2; % Raw amplitude
f=f(1:floor(N/2)); % Keep positive frequencies
amp=amp(1:floor(N/2)); % Corresponding amplitudes
figure;
plot(f,amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude/V');
title('FFT of Signal for voice 3');
hold on;
plot(freqF3,ampF3);
legend('FFT','SR760');






