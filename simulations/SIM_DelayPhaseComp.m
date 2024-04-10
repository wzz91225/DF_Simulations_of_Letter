% clc;
clear;
close all;

% Constant Definition
samp_rate = 3.2e6;  % Sampling Rate (Hz)
c = 299792458;      % Speed of light

% Drawing parameters
subplot_N = 5;      % figure 1 Number of subgraphs



% ##########################Signal##########################
% Definition of sine signal
frequency = 3.2e4;  % Frequency (Hz)
amplitude = 1;      % amplitude
sim_duration = 1;                   % Sinusoidal signal duration (s)
% Generate signal
[signal, t] = FUNC_GenerateSineSignal(frequency, amplitude, sim_duration, samp_rate);



% ##########################Gaussian noise addition##########################
% Definition of noise parameters
snr_value = 0;      % SNR (dB)
% Add noise to signal
[signal_noisy, ~] = FUNC_AddGaussianNoise(signal, snr_value);

% Draw signal
figure(1);
subplot(subplot_N, 1, 1);
% Draw original signal
plot(t, signal);
hold on;
% Draw noisy signal
plot(t, signal_noisy);
hold off;
legend('original signal', 'noisy signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
title(['signal (SNR = ' num2str(snr_value) ' dB)']);
xlim([0 100/frequency]);
ylim([-5 5]);
grid on;



% ##########################Simulation of motion delay interception##########################
% Definition of simulation related parameters
distance_relative = 20*c/frequency;     % Relative distance
relative_DoA = 66                       % Relative angle
velocity_t = 1e4;                       % Satellite motion speed
baseline_coefficient = 16;              % Interferometer baseline length coefficient (default to 2)
delta_t = c/frequency/velocity_t/baseline_coefficient;  % Comparison time interval
sampling_points_retain = samp_rate/frequency * 100;     % Comparing the number of retained sampling points

% Calculate the number of sampling point delays
alpha_sin = sin(relative_DoA * pi / 180);
alpha_cos = cos(relative_DoA * pi / 180);
dis_at = velocity_t * delta_t;
distance_A = sqrt((distance_relative * alpha_cos - dis_at/2)^2 + (distance_relative * alpha_sin)^2);
distance_B = sqrt((distance_relative * alpha_cos + dis_at/2)^2 + (distance_relative * alpha_sin)^2);
delay_A = round((distance_A / c) * samp_rate);
delay_B = round((distance_B / c + delta_t) * samp_rate);
% Intercept two sampling signals
signal_rxA = signal_noisy(delay_A : delay_A+sampling_points_retain);
signal_rxB = signal_noisy(delay_B : delay_B+sampling_points_retain);
t_rxA = t(delay_A : delay_A+sampling_points_retain);
t_rxB = t(delay_B : delay_B+sampling_points_retain);

% Draw two received intercepted signals
figure(1);
% Receiving signal interception A
subplot(subplot_N, 1, 2);
plot(t_rxA, signal_rxA);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Receiving signal interception A');
xlim([t_rxA(1) t_rxA(end)]);
ylim([-5 5]);
grid on;
% Receiving signal interception B
subplot(subplot_N, 1, 3);
plot(t_rxB, signal_rxB);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Receiving signal interception B');
xlim([t_rxB(1) t_rxB(end)]);
ylim([-5 5]);
grid on;



% ##########################Bandpass filtering##########################
% filtering
[sigA_filtered, filter_b] = FUNC_BandpassFilter(signal_rxA, frequency, samp_rate);
[sigB_filtered, ~] = FUNC_BandpassFilter(signal_rxB, frequency, samp_rate);

% View the frequency response of the filter
% figure(2)
% freqz(filter_b, 1, 1024, samp_rate);

% Draw filtered signal
figure(1);
% Receiving signal interception A
subplot(subplot_N, 1, 4);
plot(t_rxA, sigA_filtered);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Receiving signal interception A');
xlim([t_rxA(1) t_rxA(end)]);
ylim([-2 2]);
grid on;
% Receiving signal interception B
subplot(subplot_N, 1, 5);
plot(t_rxB, sigB_filtered);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Receiving signal interception B');
xlim([t_rxB(1) t_rxB(end)]);
ylim([-2 2]);
grid on;



% ##########################Time delay phase comparison direction finding##########################
distance = velocity_t * delta_t;
[~, doa_angle] = FUNC_DF2D_SignalDelayPhaseComparing(sigB_filtered, sigA_filtered, frequency, delta_t, distance, c);
doa_angle
