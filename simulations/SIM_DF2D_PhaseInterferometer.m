% clc;
clear;
close all;

% Simulation timing starts
tic;


% ##########################Simulation control##########################
is_AddNoise = 0;



% ##########################Parameter definition##########################
% Speed of Light (m/s)
c = 299792458;

% Source Frequency (Hz)
frequency = 315e6;
% Receiver signal sampling rate (Hz)
samp_rate = 315e9;

% Relative angle alpha range between signal source and receiver [0, 180)
alpha_angle = 112.5;
% Relative distance between signal source and receiverd_r (m)
d_relative = 20 * c / frequency;    % 20 times sine signal wavelength

% Definition of Gaussian noise parameters
snr_value = 30;     % SNR (dB)

% Interferometer baseline length
baseline_length = c / frequency / 2;

% Receiver sampling sine signal period number
sampling_cycles = 10;
% Total signal sampling time (s)
sampling_duration = sampling_cycles / frequency;



% ##########################Simulation space parameter calculation##########################
% Relative angle alpha related parameter
alpha_radian = deg2rad(alpha_angle);    % radian
alpha_sin = sin(alpha_radian);
alpha_cos = cos(alpha_radian);

% Spatial coordinates
[x_A, y_A] = deal(baseline_length, 0);
[x_B, y_B] = deal(-baseline_length, 0);
[x_S, y_S] = deal(d_relative * alpha_cos, d_relative * alpha_sin);

dis_AB = baseline_length;
dis_SA = sqrt((x_S - x_A)^2 + (y_S - y_A)^2);
dis_SB = sqrt((x_S - x_B)^2 + (y_S - y_B)^2);



% ##########################Simulation Time Domain Parameter Calculation##########################
% Simulation duration (s)
sim_duration = sampling_duration;
% The signal source propagates to the starting time of the first sampling at the receiver
t_start = max(dis_SA, dis_SB) / c;
% The signal source propagates to the end time of the first sampling at the receiver
t_end = max(dis_SA, dis_SB) / c + sim_duration;
% Simulation time interval (s)
sim_time_interval = 1 / samp_rate;
% Simulation time vector (s)
time_vector = t_start : sim_time_interval : t_end;



% ##########################Real time signal reception simulation##########################
% Initialize the received signal array
sig_A = zeros(1, length(time_vector));
sig_B = zeros(1, length(time_vector));

% Simulate the signal received by the receiver at each time point
for i = 1 : length(time_vector)
    
    % Calculate signal propagation time
    propagation_time_A = dis_SA / c ;
    propagation_time_B = dis_SB / c ;

    % Receiving signals
    sig_A(i) = sin(2 * pi * frequency * ...
        (time_vector(i) - propagation_time_A));
    sig_B(i) = sin(2 * pi * frequency * ...
        (time_vector(i) - propagation_time_B));
end



% ##########################Gaussian noise addition##########################
if is_AddNoise
% Add noise to signal
    sig_A_noisy = FUNC_AddGaussianNoise(sig_A, snr_value);
    sig_B_noisy = FUNC_AddGaussianNoise(sig_B, snr_value);
else
    sig_A_noisy = sig_A;
    sig_B_noisy = sig_B;
end



% ##########################Direction finding algorithm##########################
% Phase comparison
delta_phi = FUNC_ComparePhase(sig_A_noisy, sig_B_noisy);

[~, doa_angle] = FUNC_DF2D_PhaseComparing( ...
    delta_phi, frequency, baseline_length, c);



% ##########################Output results##########################
delta_phi_angle = rad2deg(delta_phi);

fprintf('Expected angle [0, 180) = %.2f°\n', alpha_angle);
fprintf('direction finding angle [0, 180) = %.2f°\n', doa_angle);
fprintf('phase difference [-180, 180] = %.2f°\n', delta_phi_angle);



% ##########################drawing##########################
% Signal drawing unified number of points
plot_points = floor(sim_duration / sim_time_interval);


% original signal and noisy signal
figure;

subplot(2, 1, 1);
plot(time_vector(1:plot_points), sig_A(1:plot_points), ...
    'DisplayName', 'Antenna A');
hold on;
plot(time_vector(1:plot_points), sig_B(1:plot_points), ...
    'DisplayName', 'Antenna B');
hold off;
legend('show');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('original signal');
xlim([time_vector(1) time_vector(plot_points)]);
ylim([-2 2]);
grid on;

subplot(2, 1, 2);
plot(time_vector(1:plot_points), sig_A_noisy(1:plot_points), ...
    'DisplayName', 'Antenna A');
hold on;
plot(time_vector(1:plot_points), sig_B_noisy(1:plot_points), ...
    'DisplayName', 'Antenna B');
hold off;
legend('show');
xlabel('Time (seconds)');
ylabel('Amplitude');
title(['noisy signal (SNR = ' num2str(snr_value) ' dB) ']);
xlim([time_vector(1) time_vector(plot_points)]);
ymax = ceil(max(max(abs(sig_A_noisy)), max(abs(sig_B_noisy))));
ylim([-ymax ymax]);
grid on;



% Simulation timing ended
toc;
