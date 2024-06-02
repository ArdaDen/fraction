close all;
clear all;
clc;

%% FRACTIONAL DELAY FINDING WITH CROSS CORRELATION CODE 

%% Variables
Fs = 1e5;  % Sampling frequency
dt = 1/Fs; % Sampling period
Time = 0.05; % Time duration PRI
dt_div = 1000; % Division of dt for fractional delay
fc = 10000; % Operating frequency
pulse_length = 5e-3; % Length of the pulse 
t = 0:dt:Time-dt; % Time scale
delay = 9.656*dt; % Fractional delay of incoming signal

%% Signal construction
window = 1*(t>=0 & t<=pulse_length); % Original pulse width
signal = sin(2*pi*fc*t).*window; % Original sinusoidal pulse signal

window_del = 1*((t-delay)>=0 & (t-delay)<=pulse_length); % Delayed pulse width
signal_del = sin(2*pi*fc*(t-delay)).*window_del; % Delayed sinusoidal pulse signal

%% Plots of the signals
figure;
plot(t,signal,"LineWidth",1,"Color","b"); 
hold on;
plot(t, signal_del,"LineWidth",1,"Color","g");
title("Original and Delayed Sinusoidal Pulse Signals");
legend("Original Signal","Delayed Signal");
xlabel("Time(s)");
ylabel("Voltage(V)");

%% Correlations
[cross_cor,t_cor] = xcorr(signal,signal_del); % Cross correlation of the original and delayed signals
[auto_cor,t_acor] = xcorr(signal); % Auto correlation of the original signal

%% Plots of the auto correlation and cross correlation functions
figure;
plot(t_cor/Fs,cross_cor,"LineWidth",1,"Color","b");
hold on;
plot(t_acor/Fs,auto_cor,"LineWidth",1,"Color","g");
title("Auto Correlation of Original Signal and Cross Correlation of Original and Delayed Signals");
legend("Cross Correlation","Auto Correlation");
xlabel("Time(s)");
ylabel("Amplitude");

%% Finding the unit delay with peak of cross correlation
[pk,i] = max(cross_cor); % Peak value and index value of the peak for cross correlation function 
unit_delay_sample = abs(t_cor(i)); % Unit sample delay
unit_delay = unit_delay_sample*dt; % Unit delay

%% Finding the unit delay with delay function
unit_sample = finddelay(signal,signal_del); % Unit sample delay
unit = unit_sample*dt; % Unit delay

%% Delaying the signal by unit delay
window = 1*((t-unit_sample*dt)>=0 & (t-unit_sample*dt)<=pulse_length);  % Original pulse width is delayed by unit delay
signal = sin(2*pi*fc*(t-unit_sample*dt)).*window; % Original signal is delayed by unit delay

%% Finding the fractional part of the delay with mini cross correlations
correlation_values = [];
for mini_delays = -dt:dt/dt_div:dt
    window_mini_del = 1*((t-unit_sample*dt-mini_delays)>=0 & (t-unit_sample*dt-mini_delays)<=pulse_length);
    signal_mini_del = sin(2*pi*fc*(t-unit_sample*dt-mini_delays)).*window_mini_del;
    cor_value = sum(signal_mini_del.*signal_del);
    correlation_values = [correlation_values;cor_value];
end
mini_delays = -dt:dt/dt_div:dt; % Mini delays scale
[pk,i] = max(correlation_values); % Finding the maximum of mini cross correlations
fractional_part = mini_delays(i); % Fractional part
total_delay = fractional_part + unit; % Total delay

%% Results
fprintf("The unit delay found by cross correlation is %0.5f microseconds.\n",unit_delay*1e6);
fprintf("The unit delay found by unit function is %0.5f microseconds.\n",unit*1e6);
fprintf("The fractional part of the delay found by mini cross correlations is %0.5f microseconds.\n",fractional_part*1e6);
fprintf("The total delay found is %0.5f microseconds.\n",total_delay*1e6);
fprintf("The real delay is  %0.5f microseconds.\n",delay*1e6);

%% Plot of the mini cross correlations
figure;
plot(mini_delays,correlation_values,"LineWidth",1,"Color","r");
title("Mini Cross Correlation of Unit Delayed Original and Delayed Signals");
xlabel("Time(s)");
ylabel("Amplitude");