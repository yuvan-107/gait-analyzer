%close all;
function avg_correlation = process_one(file_path, fs)

T = readtable(file_path);
T = renamevars(T, ["LinearAccelerationZ_m_s_2_", "Time_s_"], ["az_m_s_2", "time"]);

acc_z = T.az_m_s_2;

x = medfilt1(acc_z, 5);
Wc = 6;

[b, a] = butter(2, Wc/(fs/2));
acc_z_filt = filtfilt(b, a, x);

norm_acc = (acc_z_filt - mean(acc_z_filt)) / std(acc_z_filt);

[ACF_1, shifts] = xcorr(norm_acc, 'normalized');

p_shifts = shifts(shifts > 0)'; %taking transpose cuz its row vector
p_ACF = ACF_1(shifts > 0); %locs(1) is when shift l = 0 so take next highest. But shifts > 0 does the job for us

[pks, locs] = findpeaks(p_ACF, p_shifts, "MinPeakDistance", 120);

CNT = locs(1);

num_steps = floor(length(norm_acc) / CNT);
STEP_COUNT = reshape(norm_acc(1:num_steps*CNT), CNT, num_steps);

corr_vals = [];

for i = 2:num_steps
    temp = corrcoef(STEP_COUNT(:,i-1), STEP_COUNT(:,i));
    corr_vals(end+1) = temp(1,2);
end

avg_correlation = mean(corr_vals);

end

fs = 100;

NW_folder = "NW";
UW_folder = "UW";

NW_files = dir(fullfile(NW_folder, "*.csv"));
UW_files = dir(fullfile(UW_folder, "*.csv"));

NW_vals = [];
UW_vals = [];

for k = 1:length(NW_files)
    file_path = fullfile(NW_folder, NW_files(k).name);
    avg_corr = process_one(file_path, fs);
    NW_vals(end+1) = avg_corr;
end

for k = 1:length(UW_files)
    file_path = fullfile(UW_folder, UW_files(k).name);
    avg_corr = process_one(file_path, fs);
    UW_vals(end+1) = avg_corr;
end

threshold = (median(NW_vals) + median(UW_vals)) / 2


TP = 0; FP = 0; FN = 0; TN = 0;

for i = 1:length(NW_vals)
    if NW_vals(i) > threshold
        TP = TP + 1;
    else
        FN = FN + 1;
    end
end

for i = 1:length(UW_vals)
    if UW_vals(i) > threshold
        FP = FP + 1;
    else
        TN = TN + 1;
    end
end

precision = TP / (TP + FP)
recall = TP / (TP + FN)
F1 = 2 * (precision * recall) / (precision + recall)

% figure;
% histogram(NW_vals, 'Normalization','probability')
% hold on
% histogram(UW_vals, 'Normalization','probability')
% xline(threshold, 'r', 'LineWidth', 2)
% legend("NW","UW","Threshold")
% title("Distribution of Correlation")
% xlabel("Correlation")
% ylabel("Probability")

fprintf("\nConfusion Matrix:\n");
fprintf("[ TP = %d, FP = %d\n", TP, FP);
fprintf("  FN = %d, TN = %d ]\n\n", FN, TN);

