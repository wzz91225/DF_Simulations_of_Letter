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
% Receiver signal sampling rate (Hz)
samp_rate = 6.4e6;

% Relative distance between signal source and receiver during phase comparison d_r (m)
d_relative = 20 * c / frequency;    % 20 times sine signal wavelength
% Receiver horizontal movement speed (m/s)
v_rx = 10e3;

% Receiver phase comparison coherent accumulation sequence number
coherent_integration_number = 1;
% The receiver accumulates the number of cycles in each sequence containing a sinusoidal signal sequence compared to phase coherence
coherent_integration_cycles = 100;

% (option) Whether to use a bandpass filter (use if there is no input)
is_bandpassfilter = 0;
% % (option) Order of bandpass filter (automatically configured if there is no input or ≤ 0)
% filter_n = 200;

% ##########################Variable definition##########################
% The relative angle alpha range when comparing the signal source and receiver is [0, 180)
alpha_angle = (0:1:179);
% SNR (dB)
snr_value = (-10:10:30);

% ##########################simulation##########################
% Number of simulations
sim_num = 1000;

% Direction finding result array initialization
doa_phase_angle = ...
    zeros(length(alpha_angle), length(snr_value), sim_num);
doa_phase_angle_ave = ...
    zeros(length(alpha_angle), length(snr_value));

% Parallel loop simulation
len_alpha_angle = length(alpha_angle);
len_snr_value = length(snr_value);
% for i = 1 : len_alpha_angle
parfor i = 1 : len_alpha_angle
    for j = 1 : len_snr_value
        tmp_doa_phase_angle = 0;
        tmp_snr_value = snr_value(j);
        
        for k = 1 : sim_num
            [tmp1] = ...
                FUNC_SIM_DynamicPhaseComparingDirectionFinding( ...
                c, frequency, samp_rate, alpha_angle(i), ...
                d_relative, v_rx, tmp_snr_value, ...
                coherent_integration_number, coherent_integration_cycles, ...
                is_bandpassfilter);
            
            doa_phase_angle(i, j, k) = tmp1;

            tmp_doa_phase_angle = tmp_doa_phase_angle + tmp1;
        end
        
        doa_phase_angle_ave(i, j) = tmp_doa_phase_angle / sim_num;
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
        'SIMDATA-%s-DynamicPhaseComparingDF_Ave_SNR.mat', ...
        currentTime);
    
    % Save relevant variables to a file
    save(fileName, 'sim_num', 'sim_timing', ...
        'doa_phase_angle', ...
        'doa_phase_angle_ave', ...
        'c', 'frequency', 'samp_rate', 'alpha_angle', ...
        'd_relative', 'v_rx', 'snr_value', ...
        'coherent_integration_number', 'coherent_integration_cycles', ...
        'is_bandpassfilter');
end
