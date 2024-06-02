%% IQ DATA

%% Variables
Fs = 1e5;  % Sampling frequency
dt = 1/Fs; % Sampling period
Time = 0.05; % Time duration PRI
dt_div = 1000; % Division of dt for fractional delay
fc = 30000; % Operating frequency
pulse_length = 5e-3; % Length of the pulse 
t = 0:dt:Time-dt; % Time scale
delay = 9.656*dt; % Fractional delay of the incoming signal

%% Signal construction
window = 1*(t>=0 & t<=pulse_length); % Original pulse width
signal = sin(2*pi*fc*t).*window; % Original sinusoidal pulse signal

window_del = 1*((t-delay)>=0 & (t-delay)<=pulse_length); % Delayed pulse width
signal_del = sin(2*pi*fc*(t-delay)).*window_del; % Delayed sinusoidal pulse signal
q_data = signal_del.*cos(2*pi*fc*t); % Q Data
i_data = signal_del.*sin(2*pi*fc*t); % I Data

% Plot of Q Data
figure;
subplot(3,2,1)
plot(t,q_data,"LineWidth",1,"Color","b");
title("Q Data");
xlabel("Time(s)")
ylabel("Voltage(V)")

% Plot of I Data
subplot(3,2,2)
plot(t,i_data,"LineWidth",1,"Color","b");
title("I Data");
xlabel("Time(s)")
ylabel("Voltage(V)")

fft_q_data = fft(q_data); %FFT of Q Data
freq_scale = linspace(0,Fs,length(t)); % Frequency scale of FFT
% Plot of absolute
subplot(3,2,3)
plot(freq_scale,abs(fft_q_data)/length(t),"LineWidth",1,"Color","b");
title("FFT of Q Data")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

fft_i_data = fft(i_data); %FFT of I Data

% Plot of absolute
subplot(3,2,4)
plot(freq_scale,abs(fft_i_data)/length(t),"LineWidth",1,"Color","b");
title("FFT of I Data")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

% Plot of shifted FFT of Q Data
shifted_freq_scale = linspace(-Fs/2,Fs/2,length(t)); % Shifted frequency scale
subplot(3,2,5)
plot(shifted_freq_scale,fftshift(abs(fft_q_data)/length(t)),"LineWidth",1,"Color","b");
title("Shifted FFT of Q Data")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

% Plot of shifted FFT of I Data
subplot(3,2,6)
plot(shifted_freq_scale,fftshift(abs(fft_i_data)/length(t)),"LineWidth",1,"Color","b");
title("Shifted FFT of I Data")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

low_pass_filter = designfilt("lowpassfir","FilterOrder",16,"DesignMethod","maxflat","HalfPowerFrequency",fc*0.8,"SampleRate",Fs);
fvtool(low_pass_filter,"MagnitudeDisplay","magnitude")
q_data_filtered = filter(low_pass_filter,q_data); % Sine theta
i_data_filtered = filter(low_pass_filter,i_data); % Cosine theta 


% Plot of Q Data after filtering
figure;
subplot(3,2,1)
plot(t,q_data_filtered,"LineWidth",1,"Color","b");
title("Q Data after Lowpass Filtering")
xlabel("Time(s)")
ylabel("Voltage(V)")

% Plot of I Data after filtering
subplot(3,2,2)
plot(t,i_data_filtered,"LineWidth",1,"Color","b");
title("I Data after Lowpass Filtering")
xlabel("Time(s)")
ylabel("Voltage(V)")

% FFT of Q Data after filtering
fft_q_data_filtered = fft(q_data_filtered); % FFT

% FFT of I Data after filtering
fft_i_data_filtered = fft(i_data_filtered); % FFT

% Plot of absolute after filtering Q Data
subplot(3,2,3)
plot(freq_scale,abs(fft_q_data_filtered)/length(t),"LineWidth",1,"Color","b");
title("FFT of Q Data after Filtering")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

% Plot of absolute after filtering I Data
subplot(3,2,4)
plot(freq_scale,abs(fft_i_data_filtered)/length(t),"LineWidth",1,"Color","b");
title("FFT of I Data after Filtering")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

% Plot of shifted FFT after filtering Q Data
subplot(3,2,5)
plot(shifted_freq_scale,fftshift(abs(fft_q_data_filtered)/length(t)),"LineWidth",1,"Color","b");
title("Shifted FFT of Q Data after Filtering")
xlabel("Frequency(Hz)")
ylabel("Amplitude")

% Plot of shifted FFT after filtering I Data
subplot(3,2,6)
plot(shifted_freq_scale,fftshift(abs(fft_i_data_filtered)/length(t)),"LineWidth",1,"Color","b");
title("Shifted FFT of I Data after Filtering")
xlabel("Frequency(Hz)")
ylabel("Amplitude")  




