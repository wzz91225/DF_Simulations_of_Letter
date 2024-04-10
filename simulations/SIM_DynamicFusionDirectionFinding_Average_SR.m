% clc;
clear;
% close all;

% Initialize parallel pool (if not already started)
if isempty(gcp('nocreate'))
    % Starting parallel pooling using default settings
    parpool;
end

% Simulation timing starts
tic;

% ##########################Constant Definition##########################
% Speed of light (m/s)
c = 299792458;

% Signal source frequency (Hz)
frequency = 3.2e4;

% Relative distance between signal source and receiver during phase comparison d_r (m)
d_relative = 20 * c / frequency;    % 20 times sine signal wavelength
% Receiver horizontal movement speed (m/s)
v_rx = 10e3;

% Receiver phase comparison coherent accumulation sequence number
coherent_integration_number = 1;
% The receiver accumulates the number of cycles in each sequence containing a sinusoidal signal sequence compared to phase coherence
coherent_integration_cycles = 100;

% SNR (dB)
snr_value = -15;

% (option) Whether to use a bandpass filter (use if there is no input)
is_bandpassfilter = 0;
% % (option) Order of bandpass filter (automatically configured if there is no input or ≤ 0)
% filter_n = 200;

% ##########################Variable definition##########################
% The relative angle alpha range when comparing the signal source and receiver is [0, 180)
alpha_angle = (0:1:179);

% Receiver signal sampling rate (Hz)
samp_rate = [0.8e6 1.6e6 3.2e6 6.4e6 12.8e6 25.6e6];

% ##########################simulation##########################
% Number of simulations
sim_num = 1000;

% Direction finding result array initialization
doa_phase_angle = ...
    zeros(length(alpha_angle), length(samp_rate), sim_num);
doa_amplitude_angle = ...
    zeros(length(alpha_angle), length(samp_rate), sim_num);
doa_phase_angle_ave = ...
    zeros(length(alpha_angle), length(samp_rate));
doa_amplitude_angle_ave = ...
    zeros(length(alpha_angle), length(samp_rate));

% Parallel loop simulation
len_alpha_angle = length(alpha_angle);
len_samp_rate = length(samp_rate);
% for i = 1 : len_alpha_angle
parfor i = 1 : len_alpha_angle
    for j = 1 : len_samp_rate
        tmp_doa_phase_angle = 0;
        tmp_doa_amplitude_angle = 0;
        tmp_samp_rate = samp_rate(j);
        
        for k = 1 : sim_num
            [tmp1, tmp2] = ...
                FUNC_SIM_DynamicFusionDirectionFinding( ...
                c, frequency, tmp_samp_rate, alpha_angle(i), ...
                d_relative, v_rx, snr_value, ...
                coherent_integration_number, coherent_integration_cycles);
            
            doa_phase_angle(i, j, k) = tmp1;
            doa_amplitude_angle(i, j, k) = tmp2;

            tmp_doa_phase_angle = tmp_doa_phase_angle + tmp1;
            tmp_doa_amplitude_angle = tmp_doa_amplitude_angle + tmp2;
        end
        
        doa_phase_angle_ave(i, j) = tmp_doa_phase_angle / sim_num;
        doa_amplitude_angle_ave(i, j) = tmp_doa_amplitude_angle / sim_num;
    end
    % disp(i);
end

% Simulation timing ended
sim_timing = toc;
fprintf('Simulation duration %.6f s.\n', sim_timing);

% ##########################Save Data to File##########################
if 1
    % Get the current system time
    currentTime = datetime('now', 'TimeZone', 'local', ...
        'Format', 'yyMMdd_HHmmss');
    
    % Generate file name, including timestamp
    fileName = sprintf( ...
        'SIMDATA-%s-DynamicFusionDF_Ave_SR.mat', ...
        currentTime);
    
    % Save relevant variables to a file
    save(fileName, 'sim_num', 'sim_timing', ...
        'doa_phase_angle', 'doa_amplitude_angle', ...
        'doa_phase_angle_ave', 'doa_amplitude_angle_ave', ...
        'c', 'frequency', 'samp_rate', 'alpha_angle', ...
        'd_relative', 'v_rx', 'snr_value', ...
        'coherent_integration_number', 'coherent_integration_cycles', ...
        'is_bandpassfilter');
end
