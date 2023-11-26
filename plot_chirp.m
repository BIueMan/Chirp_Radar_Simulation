clear all
close all

% Set parameters
fs = 10000; % sampling frequency (Hz)
T = 5e-2; % duration (s)
f0 = 100; % starting frequency (Hz)
f1 = 200; % ending frequency (Hz)

% Set parameters
fs = 10000; % sampling frequency (Hz)
T = 20; % duration (s)
f0 = 10; % starting frequency (Hz)
f1 = 800; % ending frequency (Hz)

% Generate chirp signal
t = 0:1/fs:T-1/fs; % time vector
k = (f1-f0)/T; % slope of frequency modulation (Hz/s)
phi0 = 2*pi*f0*t(1); % initial phase (rad)
x = cos(2*pi*(f0*t + k/2*t.^2 + phi0)); % chirp signal

% Plot chirp signal in time domain
figure;
plot(t,x, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Chirp Signal in Time Domain');

% Compute and plot STFT
winlen = 8*128; % window length (samples)
noverlap = round(winlen/2); % overlap (samples)
nfft = winlen; % number of FFT points
figure;
spectrogram(x,winlen,noverlap,nfft,fs,'yaxis');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Chirp Signal STFT');
