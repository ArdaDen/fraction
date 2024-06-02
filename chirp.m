close all;
clear all;
clc;

Time = 5e-2; % Time length
Fs = 1e5; % Sampling frequency
dt = 1/Fs; % Sampling period
t = 0:dt:Time-dt; % Time scale
signal_length = 5e-3; % Signal length
fc = 1e5;
signal = dsp.Chirp('SweepDirection','Unidirectional','TargetFrequency',fc,'InitialFrequency',0,'TargetTime',signal_length,'SweepTime',signal_length,'SamplesPerFrame',length(t),'SampleRate',Fs);
window = 1*(t>=signal_length & t<=2*signal_length);
wave = signal()'.*window;
window1 = 1*(t>=2*signal_length & t<=3*signal_length);
wave1 = signal()'.*window1 + randn(1,length(t));
window2 = 1*(t>=0 & t<=signal_length);
wave2 = signal()'.*window2 + randn(1,length(t));
figure;
plot(t,wave);
figure;
plot(t,wave1);
figure;
plot(t,wave2);
k = signal()'.*window2;
figure;
plot(t,k)
h = wave(end:-1:1);
h = h';
x1 = conv(wave1,h);
x2 = conv(wave2,h);
tt = -2*Time:dt:2*Time;
figure;
plot(tt,x1);
figure;
plot(tt,x2);
n = fft(x1);
freq = linspace(-Fs/2,Fs/2,length(tt));
figure;
plot(freq,fftshift(abs(n)/length(tt)));