close all
clear all

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
[b, a] = butter(4, Wc/(Ws/2), 'low');
acc_z_filt = filtfilt(b, a, acc_z);

figure(2)
plot(T.time, acc_z, 'r'); hold on;
plot(T.time, acc_z_filt, 'b');
xlabel("time");
ylabel("az_m_s_2");
title("az_m_s_2 (filtered) vs time");
legend('az_m_s_2', 'az_m_s_2 after LPF');
