close all
%clear all

% T = readtable("linear_acceleration_2026-03-21_13.12.16.csv"); %test signal 1
% T = renamevars(T, ["aT_m_s_2_", "ax_m_s_2_", "ay_m_s_2_", "az_m_s_2_"], ["aT_m_s_2", "ax_m_s_2", "ay_m_s_2", "az_m_s_2"])

T = readtable("Raw Data.csv"); %test signal 2, physiobox
T = renamevars(T, ["LinearAccelerationZ_m_s_2_", "Time_s_"], ["az_m_s_2", "time"])

figure(1)
plot(T.time,T.az_m_s_2);
linkdata on;
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 vs time");
legend("show");

acc_z = T.az_m_s_2

Ws = 100;
Wc = 2.5;
[b, a] = butter(4, Wc/(Ws/2), 'low'); %TODO: tweak filter design
acc_z_filt = filtfilt(b, a, acc_z);

figure(2)
subplot(3,1,1)
plot(T.time, acc_z, 'r'); 
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 vs time")
subplot(3,1,2)
plot(T.time, acc_z_filt, 'b');
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 (filtered) vs time");

min_acc = min(acc_z_filt);
max_acc = max(acc_z_filt);

norm_acc = (acc_z_filt - min_acc) / (max_acc - min_acc); %min max normalization

subplot(3,1,3)
plot(T.time, norm_acc)
xlabel("time")
ylabel("az_m_s_2")
title("az_m_s_2 (filtered, normalized) vs time")

% norm_acc = trimdata(norm_acc,50,Side="both");
% time = trimdata(T.time, 50, Side="both");
% 
% figure(3)
% plot(time, norm_acc)
% xlabel("time")
% ylabel("az_m_s_2")
% title("az_m_s_2 (filtered, normalized, trimmed) vs time")
% legend("norm az_m_s_2")

time = T.time;

[ACF_1, shifts]= xcorr(norm_acc, 'normalized')

figure(4)
plot(shifts, ACF_1)
xlabel("Shifts")
ylabel("ACF (normalized)")
title("Autocorrelation of Normalized Acceleration")
legend("ACF")
grid on;

p_shifts = shifts(shifts > 0)'; %taking transpose cuz its row vector
p_ACF = ACF_1(shifts > 0);  %locs(1) is when shift l = 0 so take next highest. But shifts > 0 does the job for us

[pks locs] = findpeaks(p_ACF, p_shifts, "MinPeakDistance", 30); %TODO: give tweaks in 4th parameter here
CNT = locs(1)

T_full_cycle = CNT/100; %time taken to make one full gait cycle (fs = 100 Hz)

num_steps = floor(length(norm_acc)/CNT);
STEP_COUNT = reshape(norm_acc(1 : num_steps * CNT), CNT, num_steps); 
% --- ADDITION: GAIT ASYMMETRY THRESHOLD ---

% Make sure we actually found at least 2 peaks to compare
if length(pks) >= 2
    step_regularity = pks(1);   % Similarity between L and R steps
    stride_regularity = pks(2); % Similarity between L and L (or R and R) strides
    
    % Calculate Symmetry Ratio
    symmetry_ratio = step_regularity / stride_regularity;
    
    fprintf('Step Regularity (Peak 1): %.3f\n', step_regularity);
    fprintf('Stride Regularity (Peak 2): %.3f\n', stride_regularity);
    fprintf('Calculated Symmetry Ratio: %.3f\n', symmetry_ratio);
    
    % Define the threshold for natural walking
    symmetry_threshold = 0.85; 
    
    % Classification
    if symmetry_ratio >= symmetry_threshold
        disp('RESULT: Natural/Symmetrical Walking Detected.');
    else
        disp('RESULT: Unnatural/Asymmetrical Walking Detected.');
    end
else
    disp('ERROR: Could not find enough peaks to analyze symmetry. Check filter or MinPeakDistance.');
end

