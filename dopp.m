close all;
clear all;
clc;

%% Target radial velocity calculation with pulse doppler processing
%% Variables
Fs = 1e8;  % Sampling frequency
dt = 1/Fs; % Sampling period
Range_1 = 5837.463; % Range of the target 1
Range_2 = 10786.409; % Range of the target 2
min_range = 200; % Minimum range
max_range = 10000; % Maximum range
light_speed = physconst("LightSpeed"); % Light speed
pulse_length = 2*min_range/light_speed; % pulse length
pri= 2*max_range/light_speed; % PRI
prf = 1/pri; % PRF
fc = 1e9; % Operating frequency
pulse_number = 10; % Number of pulses 
t = 0:dt:pri-dt; % Time scale
rad_velocity_in_1 = 200; % Velocity of the target 1 
rad_velocity_in_2 = 120; % Velocity of the target 2 

%% For the target with radial velocity 1
%% Baseband pulse signals and complex radar matrix generation
radar_data_sin_1 = zeros(length(t),pulse_number); % Sine pulse signals
radar_data_cos_1 = zeros(length(t),pulse_number); % Cosine pulse signals

% Generation of radar matrix
for m = 1:pulse_number
    delay_m = -2*(Range_1-rad_velocity_in_1*m*pri)/light_speed; % Delay in signal 1
    window_received_m = 1*((t+delay_m)>=0 & (t+delay_m)<=pulse_length);
    radar_data_sin_1(:,m) = sin(4*pi*fc*(Range_1-m*rad_velocity_in_1*pri)/light_speed).*window_received_m; % Different phase due to velocity
    radar_data_cos_1(:,m) = cos(-4*pi*fc*(Range_1-m*rad_velocity_in_1*pri)/light_speed).*window_received_m; % Different phase due to velocity
end

radar_data_comp_1 = radar_data_cos_1 + j*radar_data_sin_1; % Complex radar matrix

% Plot of pulse integration result before matched filtering
figure;
plot(t,abs(pulsint(radar_data_comp_1,"coherent")));
title("Coherent adding of pulses before matched filtering velocity 1");
xlabel("Time(s)");
ylabel("Voltage(V)");

%% Matched filtering of pulses and matched filtered radar matrix generation
radar_data_matched_filtered_1 = zeros(2*length(t)-1,pulse_number); % Matched filtered signals
reference = 1*(t>=0 & t<=pulse_length); % Original transmitted signal in baseband

% Matched filtering
for x = 1:pulse_number
    filtered = conv(reference(end:-1:1).',radar_data_comp_1(:,x));
    radar_data_matched_filtered_1(:,x) = filtered;
end

% Plot of single pulse after matched filtering
tt = -pri+dt:dt:pri-dt;
figure;
plot(tt,abs(radar_data_matched_filtered_1(:,1)));
title("Single pulse after matched filtering velocity 1");
xlabel("Time(s)");
ylabel("Voltage(V)");

% Plot of pulse integration result after matched filtering 
figure;
plot(tt,abs(pulsint(radar_data_matched_filtered_1,"coherent")));
title("Coherent adding of pulses after matched filtering velocity 1");
xlabel("Time(s)");
ylabel("Voltage(V)");

%% Finding the peak index of pulse integration of matched filtered pulses with CFAR
guard_cell_number = 100; % Guard cell number
reference_cell_number = 600; % One side reference cell number
pfa = 1e-6;
N = 2*reference_cell_number;
alpha = N*(pfa^(-1/N)-1);
cfar_vec = [ones(1,reference_cell_number),zeros(1,guard_cell_number),ones(1,reference_cell_number)];
cfar_t = alpha*conv(abs(pulsint(radar_data_matched_filtered_1)),cfar_vec,'same')/N; % CFAR threshold
figure;
plot(tt,cfar_t)
hold on 
plot(tt,abs(pulsint(radar_data_matched_filtered_1)))
legend("1:cfar","2:sinyal");
k = abs(pulsint(radar_data_matched_filtered_1));
x = k.*(k>cfar_t);
x = x';
zer = zeros(size(x));
s = x>0;
zer(strfind([0,s(:)'],[0 1])) = 1;
idx = cumsum(zer).*s;
out = accumarray(idx(s)',x(s)',[],@(z){z'});
[p1_total,max1_total] = max(out{1});
out_pos_in_x = strfind(x,out{1});
max1_total = max1_total + out_pos_in_x;
%% Selecting the slow time row of peak index from matched filtered radar matrix
slow_time_1 = fft(radar_data_matched_filtered_1(max1_total,:));

%% FFT plots of selected slow time data
% Frequency scales
one_side_f = linspace(0,prf-prf/pulse_number,pulse_number);
two_side_f = linspace(-prf/2,prf/2-prf/pulse_number,pulse_number);

% Unshifted and shifted FFTs
figure;
tiledlayout(2,1);
ax1 = nexttile;
plot(one_side_f,abs(slow_time_1)/pulse_number);
title("Unshifted FFT of slow time from matched filtered radar matrix velocity 1");
xlabel("Frequency(Hz)");
ylabel("Amplitude");
ax2 = nexttile;
plot(two_side_f,fftshift(abs(slow_time_1)/pulse_number));
title("Shifted FFT of slow time from matched filtered radar matrix velocity 1");
xlabel("Frequency(Hz)");
ylabel("Amplitude");

%% Finding the frequency of the peak value and the conversion to speed
guard_cell_number = 5;
reference_cell_number = 20;
pfa = 1e-6;
N = 2*reference_cell_number;
alpha = N*(pfa^(-1/N)-1);
cfar_vec = [ones(1,reference_cell_number),zeros(1,guard_cell_number),ones(1,reference_cell_number)];
c_1 = alpha*conv(fftshift(abs(slow_time_1)/pulse_number),cfar_vec,'same')/N;
figure;
plot(c_1)
hold on 
plot(fftshift(abs(slow_time_1)/pulse_number))
legend("1:cfar","2:sinyal");
k_1 = fftshift(abs(slow_time_1)/pulse_number);
x_1 = k_1.*(k_1>c_1);
%x_1 = x_1';
zer_1 = zeros(size(x_1));
s_1 = x_1>0;
zer_1(strfind([0,s_1(:)'],[0 1])) = 1;
idx_1 = cumsum(zer_1).*s_1;
out_1 = accumarray(idx_1(s_1)',x_1(s_1)',[],@(z_1){z_1'});
[p1_total,max1_total] = max(out_1{1});
out_pos_in_x = strfind(x_1,out_1{1});
max1_total = max1_total + out_pos_in_x;
[p1,index_1] = max(fftshift(abs(slow_time_1)/pulse_number));
speed_1 = two_side_f(index_1)*light_speed/(fc*2);
fprintf("The speed of the first target is %.4f m/s to outside\n",speed_1)
%% Range and velocity resolution of radar
range_resolution = light_speed*pulse_length/2;
velocity_resolution = prf*light_speed/(2*length(one_side_f)*fc);

fprintf("The range resolution of the radar is %0.4f m\n",range_resolution)
fprintf("The velocity resolution of the radar is %0.4f m/s\n",velocity_resolution)

%% For the target with radial velocity 2
%% Baseband pulse signals and complex radar matrix generation
radar_data_sin_2 = zeros(length(t),pulse_number); % Sine pulse signals
radar_data_cos_2 = zeros(length(t),pulse_number); % Cosine pulse signals

% Generation of radar matrix
for m = 1:pulse_number
    delay_m = -2*(Range_1-rad_velocity_in_2*m*pri)/light_speed;
    window_received_m = 1*((t+delay_m)>=0 & (t+delay_m)<=pulse_length);
    radar_data_sin_2(:,m) = sin(4*pi*fc*(Range_1-m*rad_velocity_in_2*pri)/light_speed).*window_received_m;
    radar_data_cos_2(:,m) = cos(-4*pi*fc*(Range_1-m*rad_velocity_in_2*pri)/light_speed).*window_received_m;
end

radar_data_comp_2 = radar_data_cos_2 + j*radar_data_sin_2; % Complex radar matrix

% Plot of pulse integration result before matched filtering
figure;
plot(t,abs(pulsint(radar_data_comp_2,"coherent")));
title("Coherent adding of pulses before matched filtering velocity 2");
xlabel("Time(s)");
ylabel("Voltage(V)");

%% Matched filtering of pulses and matched filtered radar matrix generation
radar_data_matched_filtered_2 = zeros(2*length(t)-1,pulse_number); % Matched filtered signals

% Matched filtering
for x = 1:pulse_number
    filtered = conv(reference(end:-1:1).',radar_data_comp_2(:,x));
    radar_data_matched_filtered_2(:,x) = filtered;
end

% Plot of single pulse after matched filtering
figure;
plot(tt,abs(radar_data_matched_filtered_2(:,1)));
title("Single pulse after matched filtering velocity 2");
xlabel("Time(s)");
ylabel("Voltage(V)");

% Plot of pulse integration result after matched filtering 
figure;
plot(tt,abs(pulsint(radar_data_matched_filtered_2,"coherent")));
title("Coherent adding of pulses after matched filtering velocity 2");
xlabel("Time(s)");
ylabel("Voltage(V)");

%% Finding the peak index of pulse integration of matched filtered pulses
[p2_total,max2_total] = max(abs(pulsint(radar_data_matched_filtered_2)));

%% Selecting the slow time row of peak index from matched filtered radar matrix
slow_time_2 = fft(radar_data_matched_filtered_2(max2_total,:));

%% FFT plots of selected slow time data
% Unshifted and shifted FFTs
figure;
tiledlayout(2,1);
ax3 = nexttile;
plot(one_side_f,abs(slow_time_2)/pulse_number);
title("Unshifted FFT of slow time from matched filtered radar matrix velocity 2");
xlabel("Frequency(Hz)");
ylabel("Amplitude");
ax4 = nexttile;
plot(two_side_f,fftshift(abs(slow_time_2)/pulse_number));
title("Shifted FFT of slow time from matched filtered radar matrix velocity 2");
xlabel("Frequency(Hz)");
ylabel("Amplitude");

%% Finding the frequency of the peak value and the conversion to speed
[p2,index_2] = max(fftshift(abs(slow_time_2)/pulse_number));
speed_2 = two_side_f(index_2)*light_speed/(fc*2);
fprintf("The speed of the second target is %.4f m/s to outside\n",speed_2)

%% Range and velocity resolution of radar
range_resolution = light_speed*pulse_length/2;
velocity_resolution = prf*light_speed/(2*length(one_side_f)*fc);

fprintf("The range resolution of the radar is %0.4f m\n",range_resolution)
fprintf("The velocity resolution of the radar is %0.4f m/s\n",velocity_resolution)


