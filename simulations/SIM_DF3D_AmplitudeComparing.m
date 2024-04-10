% clc;
clear;
close all;

% Parameter settings
samp_rate = 32e6;         % sampling rate 32MHz
samp_duration = 0.001;    % Sampling time 0.001s
frequency = 32e3;         % frequency 32kHz
real_azimuth_angle = 30;  % Azimuth angle is 30 degrees
real_pitch_angle = 45;    % The pitch angle is 45 degrees
noise_amplitude = 0.5;    % Noise amplitude

% Generate a single frequency sine signal
t = 0:1/samp_rate:samp_duration-1/samp_rate;
sine_wave = sin(2 * pi * frequency * t);

% Generate three signals based on arrival angle
real_azimuth_radian = deg2rad(real_azimuth_angle);
real_pitch_radian = deg2rad(real_pitch_angle);
signal_ch1 = sine_wave * sin(real_pitch_radian) * cos(real_azimuth_radian);
signal_ch2 = sine_wave * sin(real_pitch_radian) * sin(real_azimuth_radian);
signal_ch3 = sine_wave * cos(real_pitch_radian);

% Output simulation setting actual angle
fprintf('Actual azimuth angle = %.2f°\n', real_azimuth_angle);
fprintf('Actual pitch angle = %.2f°\n', real_pitch_angle);

% Add Gaussian noise
signal_ch1 = signal_ch1 + noise_amplitude * randn(size(signal_ch1));
signal_ch2 = signal_ch2 + noise_amplitude * randn(size(signal_ch2));
signal_ch3 = signal_ch3 + noise_amplitude * randn(size(signal_ch3));

% Example of calculating signal-to-noise ratio (SNR), adjusted according to actual situation
snr_ch1 = 10 * log10(var(signal_ch1) / var(noise_amplitude * randn(size(signal_ch1))));
snr_ch2 = 10 * log10(var(signal_ch2) / var(noise_amplitude * randn(size(signal_ch2))));
snr_ch3 = 10 * log10(var(signal_ch3) / var(noise_amplitude * randn(size(signal_ch3))));

% Output signal-to-noise ratio
fprintf('SNR_ch1 = %.2f dB\n', snr_ch1);
fprintf('SNR_ch2 = %.2f dB\n', snr_ch2);
fprintf('SNR_ch3 = %.2f dB\n', snr_ch3);

% 3D amplitude comparison method for directional measurement
[azimuth_angle, pitch_angle] = FUNC_DF3D_AmplitudeComparing(signal_ch1, signal_ch2, signal_ch3, samp_rate);
fprintf('azimuth angle = %.2f°\n', azimuth_angle);
fprintf('pitch angle = %.2f°\n', pitch_angle);

% Draw time-domain signals
figure;
plot(t, signal_ch1, 'DisplayName', 'ch1');
hold on;
plot(t, signal_ch2, 'DisplayName', 'ch2');
plot(t, signal_ch3, 'DisplayName', 'ch3');
hold off;
legend('show');
title('Rx Signals');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
