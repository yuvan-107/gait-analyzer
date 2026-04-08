close all
%clear all

% T = readtable("linear_acceleration_2026-03-21_13.12.16.csv"); %test signal 1
% T = renamevars(T, ["aT_m_s_2_", "ax_m_s_2_", "ay_m_s_2_", "az_m_s_2_"], ["aT_m_s_2", "ax_m_s_2", "ay_m_s_2", "az_m_s_2"])

T = readtable(fullfile("NW", "Raw Data(3).csv")); %test signal 2, physiobox 
T = renamevars(T, ["LinearAccelerationZ_m_s_2_", "Time_s_"], ["az_m_s_2", "time"]);

fs = 100;

figure()
plot(T.time,T.az_m_s_2);
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 vs time");

acc_z = T.az_m_s_2;

N = length(acc_z);
f = (0:N/2) * fs / N;
X = abs(fft(acc_z));
X = X(1:N/2+1);

figure()
plot(f, X)
xlabel('Hz'); ylabel('Magnitude')
title('Frequencies of az_m_s_2')

x = medfilt1(acc_z, 5);
Wc = 3.5;
Ws = 100;
[b, a] = butter(8, Wc/(Ws/2));
acc_z_filt = filtfilt(b, a, acc_z);

figure()
plot(T.time, acc_z_filt, 'b');
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 (filtered) vs time");

min_acc = min(acc_z_filt);
max_acc = max(acc_z_filt);

%norm_acc = (acc_z_filt - min_acc) / (max_acc - min_acc); %min max normalisation
norm_acc = (acc_z_filt - mean(acc_z_filt)) / std(acc_z_filt); %Z score normalization

figure()
plot(T.time, norm_acc)
xlabel("time")
ylabel("az_m_s_2")
title("az_m_s_2 (filtered, normalized) vs time")

% norm_acc = trimdata(norm_acc,20,Side="both");
% time = trimdata(T.time, 20, Side="both");
% 
% figure()
% plot(time, norm_acc)
% xlabel("time")
% ylabel("az_m_s_2")
% title("az_m_s_2 (filtered, normalized, trimmed) vs time")
% legend("norm az_m_s_2")

time = T.time;

[ACF_1, shifts]= xcorr(norm_acc, 'normalized');

figure()
plot(shifts, ACF_1)
%xlim([0, 100]);
xlabel("Shifts")
ylabel("ACF value (normalized)")
title("Autocorrelation of Normalized Acceleration")
legend("ACF")
grid on;

p_shifts = shifts(shifts > 0)'; %taking transpose cuz its row vector
p_ACF = ACF_1(shifts > 0);  %locs(1) is when shift l = 0 so take next highest. But shifts > 0 does the job for us

[pks locs] = findpeaks(p_ACF, p_shifts, "MinPeakDistance", 120);
CNT = locs(1)

num_steps = floor((length(norm_acc) / CNT))
STEP_COUNT = reshape(norm_acc(1 : num_steps * CNT), CNT, num_steps); 

correlation_values = zeros(1, num_steps - 1);
for i = 2:num_steps
    temp_corr = (corrcoef(STEP_COUNT(:, i-1), STEP_COUNT(:, i)));
    correlation_values(i-1) = temp_corr(1, 2);
end
avg_correlation = mean(correlation_values);
fprintf('Average Left-Right Autocorrelation: %f\n', avg_correlation);

if avg_correlation > 0.50
    fprintf('Result: NATURAL WALKING\n');
else
    fprintf('Result: UNNATURAL WALKING\n');
end