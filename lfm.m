close all;
clear all;
clc;

%% This code illustrates the construction of a linear frequency modulated pulse signal. The signal is then matched filtered with and without noise and its frequency characteristics are plotted.

%% Variables
Time = 5e-2; % Time length
Fs = 1e5; % Sampling frequency
dt = 1/Fs; % Sampling period
pulse_length = 8e-3; % LFM pulse length
t = 0:dt:Time-dt; % Time scale
frequency_slope = 1000; % The slope of frequency increase
fc = 1000; % Initial frequency
snr = 0; % SNR

%% Signal construction
window = 1*(t>=0 & t<=pulse_length); % LFM signal width
signal = sin(2*pi*fc*t+2*pi*fc*frequency_slope*t.^2).*window; % Signal

%% Matched filtering
matched_filter = signal(end:-1:1).'; % Matched filter
signal_matched_filtered = conv(matched_filter,signal); % Matched filtering
tt = -Time+dt:dt:Time-dt; % Matched filtered signal time scale

%% FFT
signal_fft = fftshift(fft(signal)); % FFT of signal
freq_scale = linspace(-Fs/2,Fs/2-Fs/length(t),length(t)); % Frequency scale

%% Noise addition
signal_w_noise = awgn(signal,snr,"measured"); % Noise addition

%% Matched filtering of noisy signal
signal_w_noise_matched_filtered = conv(matched_filter,signal_w_noise); % Matched filtering of noisy signal

%% FFT of noisy signal
signal_w_noise_fft = fftshift(fft(signal_w_noise)); % FFT of noisy signal

%% Plots
figure;
tiledlayout(2,3);
ax1 = nexttile;
plot(t,signal);
title("LFM Signal")
xlabel("Time(s)")
ylabel("Voltage(V)")
ax2 = nexttile;
plot(tt,signal_matched_filtered);
title("LFM Signal after Matched Filtering")
xlabel("Time(s)")
ylabel("Voltage(V)")
ax3 = nexttile;
plot(freq_scale,abs(signal_fft)/length(t));
title("FFT of LFM Signal")
xlabel("Frequency(Hz)")
ylabel("Amplitude")
ax4 = nexttile;
plot(t,signal_w_noise);
title("Noisy LFM Signal")
xlabel("Time(s)")
ylabel("Voltage(V)")
ax5 = nexttile;
plot(tt,signal_w_noise_matched_filtered);
title("Noisy LFM Signal after Matched Filtering")
xlabel("Time(s)")
ylabel("Voltage(V)")
ax6 = nexttile;
plot(freq_scale,abs(signal_w_noise_fft)/length(t));
title("FFT of Noisy LFM Signal")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

