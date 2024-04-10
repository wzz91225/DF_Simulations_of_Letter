% clc;
clear;
close all;

% Parameter settings
samp_rate = 32e6;         % sampling rate 32MHz
samp_duration = 0.001;    % Sampling time 0.001s
frequency = 32e3;         % frequency 32kHz
real_azimuth_angle = 67;  % Azimuth angle (degree)
noise_amplitude = 0.5;    % Noise amplitude

% Generate a single frequency sine signal
t = 0:1/samp_rate:samp_duration-1/samp_rate;
sine_wave = sin(2 * pi * frequency * t);

% Generate two signals based on arrival angle
real_azimuth_radian = deg2rad(real_azimuth_angle);
signal_ch1 = sine_wave * cos(real_azimuth_radian); % X-axis signal
signal_ch2 = sine_wave * sin(real_azimuth_radian); % Y-axis signal

% Output simulation setting actual angle
fprintf('actual angle = %.2f°\n', real_azimuth_angle);

% Add Gaussian noise
signal_ch1 = signal_ch1 + noise_amplitude * randn(size(signal_ch1));
signal_ch2 = signal_ch2 + noise_amplitude * randn(size(signal_ch2));

% azimuth angle
[~, azimuth_angle] = FUNC_DF2D_AmplitudeComparing( ...
    signal_ch1, signal_ch2, samp_rate);

fprintf('azimuth angle = %.2f°\n', azimuth_angle);

% Draw time-domain signals
figure;
plot(t, signal_ch1, 'DisplayName', 'ch1 - X axis');
hold on;
plot(t, signal_ch2, 'DisplayName', 'ch2 - Y axis');
hold off;
legend('show');
title('Rx Signals - 2D Direction Finding');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
