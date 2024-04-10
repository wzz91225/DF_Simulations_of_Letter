% clc;
clear;
close all;

% Simulation timing starts
tic;



% ##########################Simulation process control##########################
% Whether to output results
is_fprintf = 1;
% Whether to draw
is_figure = 1;

% whether to use a filter
is_filter = 0;
% Whether to use coherent accumulators
is_coherent_integration = 1;



% ##########################Definition of Main Simulation Parameters##########################
% Speed of light (m/s)
c = 299792458;

% Signal source frequency (Hz)
frequency = 3.2e4;
% Receiver signal sampling rate (Hz)
samp_rate = 6.4e6;

% The relative angle alpha range when comparing the signal source and receiver is [0, 180)
alpha_angle = 34;
% Relative distance between signal source and receiver during phase comparison d_r (m)
d_relative = 20 * c / frequency;    % 20 times sine signal wavelength
% Receiver horizontal movement speed (m/s)
v_rx = 10e3;

% SNR (dB)
snr_value = 0;

% Receiver phase comparison coherent integration sequence number
coherent_integration_number = 10;
% The receiver integrates the number of cycles in each sequence containing a sinusoidal signal sequence compared to phase coherence
coherent_integration_cycles = 10;



% ##########################Signal Source & Receiver Parameter Calculation##########################
% Total cycle number of single sampling sinusoidal signal in receiver phase comparison
single_sampling_cycles = coherent_integration_number * ...
    coherent_integration_cycles;
% receiver phase single-sampling duration (s)
single_sampling_duration = single_sampling_cycles / frequency;
% Receiver phase comparison distance traveled by a single sample (s)
single_sampling_distance = single_sampling_duration / v_rx;

% Receiver phase comparison interval distance (m)
sampling_interval = c / frequency / 2;  % Half wavelength
% Receiver phase comparison interval time (s)
delta_t = sampling_interval / v_rx;

% Total signal sampling time (s)
sampling_duration = delta_t + single_sampling_duration;
% The distance corresponding to the total duration of signal sampling (m)
sampling_distance = sampling_duration / v_rx;

% Relative angle alpha related parameters
alpha_radian = deg2rad(alpha_angle);    % radian
alpha_sin = sin(alpha_radian);
alpha_cos = cos(alpha_radian);

% The distance d_v of the signal source projection in the direction of receiver movement  (m)
d_vertical = d_relative * alpha_sin;
% Relative distance d'_r from signal source when receiver phase comparison first sampling is completed (m)
d_prime_relative = sqrt(d_vertical^2 + ...
    (sampling_interval + d_relative * alpha_cos)^2);

% Receiver phase comparison relative distance from signal source d_start at the start of final sampling (m)
d_start = sqrt(d_vertical^2 + ...
    (sampling_interval + d_relative * alpha_cos)^2);
% Receiver phase comparison relative distance d'_start from signal source at the start of first sampling (m)
d_prime_start = sqrt(d_vertical^2 + ...
    (sampling_distance + d_relative * alpha_cos)^2);

% Receiver initial spatial coordinates
[x_rx, y_rx] = deal(0, 0);
% Signal source coordinates
[x_s, y_s] = deal(d_relative * alpha_cos + sampling_distance, ...
    d_vertical);



% ##########################Simulation Time Domain Parameter Calculation##########################
% Simulation duration (s)
sim_duration = sampling_duration;
% The signal source propagates to the starting time of the first sampling at the receiver
t_prime_start = d_prime_start / c;
% Simulation time interval (s)
sim_time_interval = 1 / samp_rate;
% Simulation time vector (s)
time_vector = t_prime_start : sim_time_interval : ...
    t_prime_start+sim_duration;

% Phase comparison single sampling points
single_sampling_points = round( ...
    single_sampling_duration / sim_time_interval);
% points of comparison intervals
interval_points =  round(delta_t / sim_time_interval);
% Coherent integration signal sampling points
coherent_integration_points = single_sampling_points / ...
    coherent_integration_number;



% ##########################Calculation of spatial grid parameters##########################
% Maximum sampling interval distance (m)
d_max = v_rx / samp_rate;
% Simulation grid length (m)
grid_width = sampling_distance;
% Number of simulation grid points
num_points_width = ceil(grid_width / d_max);



% ##########################Real time signal reception simulation##########################
% Initialize the received signal array
sig_rx = zeros(1, length(time_vector));
sig_rx_ch1 = zeros(1, length(time_vector));
sig_rx_ch2 = zeros(1, length(time_vector));

% Simulate the signal received by the receiver at each time point
for i = 1 : length(time_vector)
    % Calculate the current position of the receiver (moving along the X-axis starting from t=0)
    x_rx = x_rx + v_rx * sim_time_interval;
    y_rx = 0;
    
    % Calculate the relative distance between the receiver and the signal source
    distance = sqrt((x_rx - x_s)^2 + (y_s - y_rx)^2);
    % Calculate signal propagation time
    propagation_time = distance / c ;
    
    % Receiving signals
    sig_rx(i) = sin(2 * pi * frequency * ...
        (time_vector(i) - propagation_time));
    
    % Dual channel orthogonal antenna receiving signal
    sig_rx_ch1(i) = sig_rx(i) * abs(alpha_cos);
    sig_rx_ch2(i) = sig_rx(i) * alpha_sin;
end



% ##########################Gaussian noise addition##########################
% Add noise to signal
sig_rx_ch1_noisy = FUNC_AddGaussianNoise(sig_rx_ch1, snr_value);
sig_rx_ch2_noisy = FUNC_AddGaussianNoise(sig_rx_ch2, snr_value);



% ##########################Direction finding signal interception##########################
% The time vectors corresponding to signals A and B respectively
tv_sigA = time_vector(1 : single_sampling_points);
tv_sigB = time_vector((1 + interval_points) : ...
    (interval_points + single_sampling_points));

% Signal interception index
idx_A_head = 1;
idx_A_tail = single_sampling_points;
idx_B_head = 1 + interval_points;
idx_B_tail = interval_points + single_sampling_points;

% Signal interception
sigA_ch1 = sig_rx_ch1_noisy(idx_A_head : idx_A_tail);
sigA_ch2 = sig_rx_ch2_noisy(idx_A_head : idx_A_tail);

sigB_ch1 = sig_rx_ch1_noisy(idx_B_head : idx_B_tail);
sigB_ch2 = sig_rx_ch2_noisy(idx_B_head : idx_B_tail);



% ##########################Frequency detection##########################
% Dual channel signal addition
sigA_sum = sigA_ch1 + sigA_ch2;
sigB_sum = sigB_ch1 + sigB_ch2;

% Calculate the signal power spectrum and its corresponding frequency vector
[fv_sigA, pspectrum_sigA] = FUNC_TransForm2PowerSpectrum( ...
    sigA_sum, samp_rate);
[fv_sigB, pspectrum_sigB] = FUNC_TransForm2PowerSpectrum( ...
    sigB_sum, samp_rate);

% Search for power spectrum peaks and their corresponding frequency points
[freq_sigA, ppower_sigA] = FUNC_FindMaxPeak(fv_sigA, pspectrum_sigA);
[freq_sigB, ppower_sigB] = FUNC_FindMaxPeak(fv_sigB, pspectrum_sigB);
if isnan(freq_sigA)
    freq_sigA = freq_sigB;
    ppower_sigA = 0;
end
if isnan(freq_sigB)
    freq_sigB = freq_sigA;
    ppower_sigB = 0;
end



% ##########################Bandpass filtering##########################
if is_filter
    % filtering
    [sigA_ch1_filtered, filter_b] = FUNC_BandpassFilter( ...
        sigA_ch1, frequency, samp_rate);
    [sigA_ch2_filtered, ~] = FUNC_BandpassFilter( ...
        sigA_ch2, frequency, samp_rate);
    
    [sigB_ch1_filtered, ~] = FUNC_BandpassFilter( ...
        sigB_ch1, frequency, samp_rate);
    [sigB_ch2_filtered, ~] = FUNC_BandpassFilter( ...
        sigB_ch2, frequency, samp_rate);
else
    sigA_ch1_filtered = sigA_ch1;
    sigA_ch2_filtered = sigA_ch2;
    sigB_ch1_filtered = sigB_ch1;
    sigB_ch2_filtered = sigB_ch2;
end



% ##########################Coherent integration##########################
if is_coherent_integration
    % The time vector corresponding to the coherent integration signal
    tv_sigA_integration = tv_sigA(end-coherent_integration_points+1 : end);
    tv_sigB_integration = tv_sigB(end-coherent_integration_points+1 : end);
    
    % Coherent integration
    sigA_ch1_integration = FUNC_SignalCoherentIntegration( ...
        sigA_ch1_filtered, coherent_integration_points, coherent_integration_number);
    
    sigA_ch2_integration = FUNC_SignalCoherentIntegration( ...
        sigA_ch2_filtered, coherent_integration_points, coherent_integration_number);
    
    sigB_ch1_integration = FUNC_SignalCoherentIntegration( ...
        sigB_ch1_filtered, coherent_integration_points, coherent_integration_number);
    
    sigB_ch2_integration = FUNC_SignalCoherentIntegration( ...
        sigB_ch2_filtered, coherent_integration_points, coherent_integration_number);
else
    tv_sigA_integration = tv_sigA;
    tv_sigB_integration = tv_sigB;

    sigA_ch1_integration = sigA_ch1_filtered;
    sigA_ch2_integration = sigA_ch2_filtered;
    sigB_ch1_integration = sigB_ch1_filtered;
    sigB_ch2_integration = sigB_ch2_filtered;
end


% ##########################Direction finding algorithm##########################
% Dynamic phase comparison direction finding algorithm
sigA_integration_sum = sigA_ch1_integration + sigA_ch2_integration;
sigB_integration_sum = sigB_ch1_integration + sigB_ch2_integration;
[~, doa_phase_angle] = FUNC_DF2D_SignalDelayPhaseComparing( ...
    sigB_integration_sum, sigA_integration_sum, frequency, ...
    delta_t, sampling_interval, c);

% Amplitude comparison direction finding algorithm
[~, doa_amplitude_angle] = FUNC_DF2D_AmplitudeComparing( ...
    sigB_ch1_integration, sigB_ch2_integration, samp_rate);

% Fusion direction finding
doa_fusion_angle = FUNC_DF2D_DirectionFindingFusionModel( ...
    doa_amplitude_angle, doa_phase_angle);


% ##########################Output results##########################
if is_fprintf
    fprintf('Frequency detection A = %.2fHz\n', freq_sigA);
    fprintf('Frequency detection B = %.2fHz\n', freq_sigB);
    fprintf('Actual angle [0, 180) = %.2f°\n', alpha_angle);
    fprintf('Angle of amplitude comparison algorithm [0, 90] = %.2f°\n', doa_amplitude_angle);
    fprintf('Angle of dynamic phase comparison algorithm [0, 180) = %.2f°\n', doa_phase_angle);
    fprintf('Angle of dynamic fusion method [0, 180) = %.2f°\n', doa_fusion_angle);
end



% ##########################drawing##########################
if is_figure
    % close all;
    % Signal drawing unified number of points
    plot_points = coherent_integration_points;


    % Original signal and noisy signal
    figure;

    subplot(2, 1, 1);
    plot(time_vector(1:plot_points), sig_rx_ch1(1:plot_points), ...
        'DisplayName', 'X-axis channel 1');
    hold on;
    plot(time_vector(1:plot_points), sig_rx_ch2(1:plot_points), ...
        'DisplayName', 'Y-axis channel 2');
    plot(time_vector(1:plot_points), sig_rx(1:plot_points), ...
        'DisplayName', 'original signal');
    hold off;
    legend('show');
    xlabel('time (s)');
    ylabel('amplitude');
    title('original signal');
    xlim([time_vector(1) time_vector(plot_points)]);
    ylim([-2 2]);
    grid on;

    subplot(2, 1, 2);
    plot(time_vector(1:plot_points), sig_rx_ch1_noisy(1:plot_points), ...
        'DisplayName', 'X-axis channel 1');
    hold on;
    plot(time_vector(1:plot_points), sig_rx_ch2_noisy(1:plot_points), ...
        'DisplayName', 'Y-axis channel 2');
    hold off;
    legend('show');
    xlabel('time (s)');
    ylabel('amplitude');
    title(['Gaussian noise signal from dual channel antenna (SNR = ' num2str(snr_value) ' dB) ']);
    xlim([time_vector(1) time_vector(plot_points)]);
    ymax = ceil(max(max(abs(sig_rx_ch1_noisy)), max(abs(sig_rx_ch2_noisy))));
    ylim([-ymax ymax]);
    grid on;


    % Dual channel interception signals A and B
    figure;

    subplot(2, 1, 1);
    plot(tv_sigA(1:plot_points), sigA_ch1(1:plot_points), ...
        'DisplayName', 'X-axis channel 1');
    hold on;
    plot(tv_sigA(1:plot_points), sigA_ch2(1:plot_points), ...
        'DisplayName', 'Y-axis channel 2');
    hold off;
    legend('show');
    xlabel('time (s)');
    ylabel('amplitude');
    title('Dual channel antenna intercepts received signal A');
    xlim([tv_sigA(1) tv_sigA(plot_points)]);
    ymax = ceil(max(max(abs(sigA_ch1)), max(abs(sigA_ch2))));
    ylim([-ymax ymax]);
    grid on;

    subplot(2, 1, 2);
    plot(tv_sigB(1:plot_points), sigB_ch1(1:plot_points), ...
        'DisplayName', 'X-axis channel 1');
    hold on;
    plot(tv_sigB(1:plot_points), sigB_ch2(1:plot_points), ...
        'DisplayName', 'Y-axis channel 2');
    hold off;
    legend('show');
    xlabel('time (s)');
    ylabel('amplitude');
    title('Dual channel antenna intercepts received signal B');
    xlim([tv_sigB(1) tv_sigB(plot_points)]);
    ymax = ceil(max(max(abs(sigB_ch1)), max(abs(sigB_ch2))));
    ylim([-ymax ymax]);
    grid on;

    
    % Dual channel interception of signal A and B spectrum
    figure;

    subplot(2, 1, 1);
    plot(fv_sigA, pspectrum_sigA);
    xlabel('frequency (Hz)');
    ylabel('amplitude');
    title('Dual channel antenna intercepts the received signal A spectrum');
    grid on;

    subplot(2, 1, 2);
    plot(fv_sigB, pspectrum_sigB);
    xlabel('frequency (Hz)');
    ylabel('amplitude');
    title('Dual channel antenna intercepts the received signal B spectrum');
    grid on;

    
    if is_filter
        % Filter frequency response and phase response
        figure;
        freqz(filter_b, 1, 1024, samp_rate);
    
    
        % 带通滤波信号A和B
        figure;
    
        subplot(2, 1, 1);
        plot(tv_sigA(1:plot_points), sigA_ch1_filtered(1:plot_points), ...
            'DisplayName', 'X-axis channel 1');
        hold on;
        plot(tv_sigA(1:plot_points), sigA_ch2_filtered(1:plot_points), ...
            'DisplayName', 'Y-axis channel 2');
        hold off;
        legend('show');
        xlabel('time (s)');
        ylabel('amplitude');
        title('Bandpass filtered dual channel signal A');
        xlim([tv_sigA(1) tv_sigA(plot_points)]);
        ylim([-2 2]);
        grid on;
    
        subplot(2, 1, 2);
        plot(tv_sigB(1:plot_points), sigB_ch1_filtered(1:plot_points), ...
            'DisplayName', 'X-axis channel 1');
        hold on;
        plot(tv_sigB(1:plot_points), sigB_ch2_filtered(1:plot_points), ...
            'DisplayName', 'Y-axis channel 2');
        hold off;
        legend('show');
        xlabel('time (s)');
        ylabel('amplitude');
        title('Bandpass filtered dual channel signal B');
        xlim([tv_sigB(1) tv_sigB(plot_points)]);
        ylim([-2 2]);
        grid on;
    end


    if is_coherent_integration
        % Coherent integration signals A and B
        figure;
    
        subplot(2, 1, 1);
        plot(tv_sigA_integration, sigA_ch1_integration, ...
            'DisplayName', 'X-axis channel 1');
        hold on;
        plot(tv_sigA_integration, sigA_ch2_integration, ...
            'DisplayName', 'Y-axis channel 2');
        hold off;
        legend('show');
        xlabel('time (s)');
        ylabel('amplitude');
        title('Coherent integration dual channel signal A');
        xlim([tv_sigA_integration(1) tv_sigA_integration(end)]);
        ylim([-2 2]);
        grid on;
    
        subplot(2, 1, 2);
        plot(tv_sigB_integration, sigB_ch1_integration, ...
            'DisplayName', 'X-axis channel 1');
        hold on;
        plot(tv_sigB_integration, sigB_ch2_integration, ...
            'DisplayName', 'Y-axis channel 2');
        hold off;
        legend('show');
        xlabel('time (s)');
        ylabel('amplitude');
        title('Coherent integration dual channel signal B');
        xlim([tv_sigB_integration(1) tv_sigB_integration(end)]);
        ylim([-2 2]);
        grid on;
    end
end

% Simulation timing ended
toc;
