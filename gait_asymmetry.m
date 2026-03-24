close all
%clear all

T = readtable("linear_acceleration_2026-03-21_13.12.16.csv");
T = renamevars(T, ["aT_m_s_2_", "ax_m_s_2_", "ay_m_s_2_", "az_m_s_2_"], ["aT_m_s_2", "ax_m_s_2", "ay_m_s_2", "az_m_s_2"])

figure(1)
plot(T.time,T.az_m_s_2,XDataSource = 'T.time',YDataSource = 'T.az_m_s_2');
linkdata on;
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 vs time");
legend("show");

acc_z = T.az_m_s_2

Ws = 100;
Wc = 3;
[b, a] = butter(4, Wc/(Ws/2), 'low'); %TODO: tweak filter design
acc_z_filt = filtfilt(b, a, acc_z);

figure(2)
plot(T.time, acc_z, 'r'); hold on;
plot(T.time, acc_z_filt, 'b');
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 (filtered) vs time");
legend('az_m_s_2', 'az_m_s_2 after LPF');

min_acc = min(acc_z_filt);
max_acc = max(acc_z_filt);

norm_acc = (acc_z_filt - min_acc) / (max_acc - min_acc); %min max normalization

figure(3)
plot(T.time, norm_acc)
xlabel("time")
ylabel("az_m_s_2")
title("az_m_s_2 (filtered, normalized) vs time")
legend("norm az_m_s_2")

% norm_acc = trimdata(norm_acc,50,Side="both");
% time = trimdata(T.time, 50, Side="both");
% 
% figure(4)
% plot(time, norm_acc)
% xlabel("time")
% ylabel("az_m_s_2")
% title("az_m_s_2 (filtered, normalized, trimmed) vs time")
% legend("norm az_m_s_2")

time = T.time;

[ACF_1, shifts]= xcorr(norm_acc, 'normalized')

figure(5)
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

t_full_cycle = CNT/100; %time taken to make one full gait cycle (fs = 100 Hz)



