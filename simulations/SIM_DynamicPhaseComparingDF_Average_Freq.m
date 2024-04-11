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

% ##########################Parameter definition##########################
% The relative angle alpha range when comparing the signal source and receiver is [0, 180)
alpha_angle = (0:1:179);

% Speed of light (m/s)
c = 299792458;


% Signal source frequency (Hz)
frequency = [3e3 3e4 3e5 3e6 3e7];
% Receiver signal sampling rate (Hz)
samp_rate = frequency * 200;
% Relative distance between signal source and receiver during phase comparison d_r (m)
d_relative = 20 * c ./ frequency;    % 20 times sine signal wavelength


% Receiver horizontal movement speed (m/s)
v_rx = 10e3;

% SNR (dB)
snr_value = -15;

% Receiver phase comparison coherent accumulation sequence number
coherent_integration_number = 1;
% The receiver accumulates the number of cycles in each sequence containing a sinusoidal signal sequence compared to phase coherence
coherent_integration_cycles = 100;

% (option) 是否使用带通滤波器（无输入则使用）
is_bandpassfilter = 0;
% % (option) 带通滤波器阶数（无输入或≤0则自动配置）
% filter_n = 200;

% ##########################simulation##########################
% Number of simulations
sim_num = 1000;

% Direction finding result array initialization
doa_phase_angle = ...
    zeros(length(alpha_angle), length(frequency), sim_num);
doa_phase_angle_ave = ...
    zeros(length(alpha_angle), length(frequency));

% Parallel loop simulation
len_alpha_angle = length(alpha_angle);
len_frequency = length(frequency);
% for i = 1 : len_alpha_angle
parfor i = 1 : len_alpha_angle
    for j = 1 : len_frequency
        tmp_doa_phase_angle = 0;
        tmp_frequency = frequency(j);
        tmp_samp_rate = samp_rate(j);
        tmp_d_relative = d_relative(j);
        
        
        for k = 1 : sim_num
            [tmp1] = ...
                FUNC_SIM_DynamicPhaseComparingDirectionFinding( ...
                c, tmp_frequency, tmp_samp_rate, alpha_angle(i), ...
                tmp_d_relative, v_rx, snr_value, ...
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
        'SIMDATA-%s-DynamicPhaseComparingDF_Ave_Freq.mat', ...
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
