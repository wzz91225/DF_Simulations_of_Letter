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

average_num = 1;

snr_value = 30;
sampling_periods = 10;

baseline_coefficient = 2;
frequency = 640e3;
samp_rate = 6.4e6;

% The relative angle alpha range when comparing the signal source and receiver is [0, 180)
relative_DoA = (0:1:179);
doa_angle = zeros(1, length(relative_DoA));

i = 0;
% for relative_DoA = 0:1:179
parfor i = 1 : length(relative_DoA)
    sum = 0;
    for j = 1:average_num
        sum = sum + FUNC_SIM_DelayPhaceComp( ...
            relative_DoA(i), snr_value, sampling_periods, ...
            baseline_coefficient, samp_rate, frequency);
    end
    doa_angle(i) = abs(sum / average_num);
end


% Simulation timing ended
sim_timing = toc;
fprintf('仿真历时 %.6f 秒。\n', sim_timing);


figure;
plot(relative_DoA, doa_angle);
hold on;
plot((0:1:179));
hold off;
axis([0 180 0 180]);
grid on;
