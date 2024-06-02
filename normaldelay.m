close all;
clear all;
clc;
Fs = 1e5;
dt = 1/Fs;
Time = 0.05;
fc = 1000;
t = -Time:dt:Time;
delay = 10*dt;
window = 1*(t>=0 & t<=5e-3);
signal = sin(2*pi*fc*t).*window;
%signal = signal';
plot(t,signal);
window2 = 1*((t-delay)>=0 & (t-delay)<=5e-3);
signal2 = sin(2*pi*fc*(t-delay)).*window2;
%signal2 = signal2';
figure;
plot(t,signal2);
hold on;
plot(t,signal);



close all;
clear all;
clc;
Fs = 1e4;
dt = 1/Fs;
Time = 1e-2;
t = 0:dt:Time;
delay = 5*dt;
window = 1*(t>=0 & t<=5e-3);
signal = sin(2*pi*700*t).*window;
signal = signal';
plot(t,signal);
FD = 0.56;
N = 10;
h = designFracDelayFIR(FD,N);
fdfir = dsp.FIRFilter(h);
y = fdfir(signal);
plot(t,signal);
hold on;
plot(t,y);
window = 1*((t-delay)>=0 & (t-delay)<=5e-3);
signal = sin(2*pi*700*(t-delay)).*window;
signal = signal';
figure;
plot(t,signal)


close all;
clear all;
clc;
Fs = 1e4;
dt = 1/Fs;
Time = 0.05;
fc = 1000;
t = -Time:dt:Time;
delay = 5.4*dt;
window = 1*(t>=0 & t<=3e-3);
signal = sin(2*pi*fc*t).*window;
signal = signal';
plot(t,signal);
delayed = delayseq(signal,6.7*dt,Fs);
finddelay(signal,delayed)
figure;
plot(t,signal);
hold on;
plot(t,delayed);
[cor,x] = xcorr(signal,delayed);
[cor2,xp] = xcorr(signal);
figure;
plot(x/Fs,cor);
hold on;
plot(x/Fs,cor2);
[pk,i] = max(cor);
pk
x(i)/Fs