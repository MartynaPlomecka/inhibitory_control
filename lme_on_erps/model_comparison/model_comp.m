addpath \\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_on_erps\stimulus_locked
load('all_subjects_ERP_lme_table_separated_direction.mat')

all_subjects = all_subject_tables.SubjectID;
all_ag_rep = all_subject_tables.AgeGroup;    % AgeGroup: young/old
all_cond = all_subject_tables.Condition;     % Condition: pro/anti
all_dir = all_subject_tables.Direction;      % Direction: left/right
all_mean_rt = all_subject_tables.MeanRT;
all_eeg_data = all_subject_tables.EEGData;   % Should be [channels x time points] matrices

n_samples = length(all_subjects);
[channels, times] = size(all_eeg_data{1});

time_points = linspace(-800, 798, times);

%  [n_samples x channels x times]
all_eeg = nan(n_samples, channels, times);
for i = 1:n_samples
    all_eeg(i, :, :) = all_eeg_data{i};
end

%% Data normalization and time window selection
time_window = [-100, 500]; % ms
window_indices = find(time_points >= time_window(1) & time_points <= time_window(2));

all_eeg = all_eeg(:, :, window_indices);
times = length(window_indices);

time_points = time_points(window_indices);

all_eeg_mean = mean(all_eeg(:));
all_eeg_std = std(all_eeg(:));

all_eeg_normalized = (all_eeg - all_eeg_mean) / all_eeg_std;

mean_rt_mean = mean(all_mean_rt);
mean_rt_std = std(all_mean_rt);
all_mean_rt_normalized = (all_mean_rt - mean_rt_mean) / mean_rt_std;

%% Mixed models on each time-channel combination

AgeGroup = categorical(all_ag_rep, [0, 1], {'young', 'old'});
AgeGroup = reordercats(AgeGroup, {'young', 'old'});

Condition = categorical(all_cond, {'pro', 'anti'}, 'Ordinal', false);
Condition = reordercats(Condition, {'pro', 'anti'});

Direction = categorical(all_dir, {'left', 'right'}, 'Ordinal', false);
Direction = reordercats(Direction, {'left', 'right'});

SubjectID = categorical(all_subjects);

%% Initialize arrays to store p-values and AIC/BIC values
p_values_1_vs_2_temp = cell(1, channels * times);
p_values_2_vs_3_temp = cell(1, channels * times);
p_values_3_vs_4_temp = cell(1, channels * times);
p_values_1_vs_4_temp = cell(1, channels * times);
p_values_1_vs_3_temp = cell(1, channels * times); % New comparison
p_values_2_vs_4_temp = cell(1, channels * times); % New comparison

AIC_values_temp = cell(1, channels * times);
BIC_values_temp = cell(1, channels * times);

parfor idx = 1:channels * times
    [ch, t] = ind2sub([channels, times], idx);

    eeg_data_for_point = squeeze(all_eeg_normalized(:, ch, t));

    eeg_tbl = table(eeg_data_for_point, AgeGroup, SubjectID, ...
        Condition, Direction, all_mean_rt_normalized, ...
        'VariableNames', {'EEG', 'AgeGroup', 'SubjectID', 'Condition', 'MeanRT', 'Direction'});

    AIC_temp = nan(4, 1);
    BIC_temp = nan(4, 1);

    lme_model_1 = fitlme(eeg_tbl, 'EEG ~ AgeGroup + Condition + MeanRT + Direction + (1|SubjectID)');
    lme_model_2 = fitlme(eeg_tbl, 'EEG ~ AgeGroup * Condition + MeanRT + Direction + (1|SubjectID)');
    lme_model_3 = fitlme(eeg_tbl, 'EEG ~ AgeGroup * Condition * MeanRT + Direction + (1|SubjectID)');
    lme_model_4 = fitlme(eeg_tbl, 'EEG ~ AgeGroup + Condition * MeanRT + Direction + (1|SubjectID)');

    AIC_temp(1) = lme_model_1.ModelCriterion.AIC;
    BIC_temp(1) = lme_model_1.ModelCriterion.BIC;

    AIC_temp(2) = lme_model_2.ModelCriterion.AIC;
    BIC_temp(2) = lme_model_2.ModelCriterion.BIC;

    AIC_temp(3) = lme_model_3.ModelCriterion.AIC;
    BIC_temp(3) = lme_model_3.ModelCriterion.BIC;

    AIC_temp(4) = lme_model_4.ModelCriterion.AIC;
    BIC_temp(4) = lme_model_4.ModelCriterion.BIC;

    % Comparing models using likelihood ratio tests and store p-values temporarily
    results_1_vs_2 = compare(lme_model_1, lme_model_2);
    results_2_vs_3 = compare(lme_model_2, lme_model_3);
    results_3_vs_4 = compare(lme_model_3, lme_model_4);
    results_1_vs_4 = compare(lme_model_1, lme_model_4);
    results_1_vs_3 = compare(lme_model_1, lme_model_3);
    results_2_vs_4 = compare(lme_model_2, lme_model_4);

    p_values_1_vs_2_temp{idx} = results_1_vs_2.pValue(2);
    p_values_2_vs_3_temp{idx} = results_2_vs_3.pValue(2);
    p_values_3_vs_4_temp{idx} = results_3_vs_4.pValue(2);
    p_values_1_vs_4_temp{idx} = results_1_vs_4.pValue(2);
    p_values_1_vs_3_temp{idx} = results_1_vs_3.pValue(2);
    p_values_2_vs_4_temp{idx} = results_2_vs_4.pValue(2);

    % Store AIC and BIC for later analysis
    AIC_values_temp{idx} = AIC_temp;
    BIC_values_temp{idx} = BIC_temp;
end

%% Reconstruct the p-value matrices after the parfor loop
p_values_1_vs_2 = nan(channels, times);
p_values_2_vs_3 = nan(channels, times);
p_values_3_vs_4 = nan(channels, times);
p_values_1_vs_4 = nan(channels, times);
p_values_1_vs_3 = nan(channels, times); % New comparison
p_values_2_vs_4 = nan(channels, times); % New comparison

for idx = 1:channels * times
    [ch, t] = ind2sub([channels, times], idx);
    p_values_1_vs_2(ch, t) = p_values_1_vs_2_temp{idx};
    p_values_2_vs_3(ch, t) = p_values_2_vs_3_temp{idx};
    p_values_3_vs_4(ch, t) = p_values_3_vs_4_temp{idx};
    p_values_1_vs_4(ch, t) = p_values_1_vs_4_temp{idx};
    p_values_1_vs_3(ch, t) = p_values_1_vs_3_temp{idx};
    p_values_2_vs_4(ch, t) = p_values_2_vs_4_temp{idx};
end

%% Reconstruct the AIC and BIC matrices from the temporary cell arrays
AIC_values_model_1 = nan(channels, times);
AIC_values_model_2 = nan(channels, times);
AIC_values_model_3 = nan(channels, times);
AIC_values_model_4 = nan(channels, times);

BIC_values_model_1 = nan(channels, times);
BIC_values_model_2 = nan(channels, times);
BIC_values_model_3 = nan(channels, times);
BIC_values_model_4 = nan(channels, times);

for idx = 1:channels * times
    [ch, t] = ind2sub([channels, times], idx);
    AIC_values_model_1(ch, t) = AIC_values_temp{idx}(1);
    AIC_values_model_2(ch, t) = AIC_values_temp{idx}(2);
    AIC_values_model_3(ch, t) = AIC_values_temp{idx}(3);
    AIC_values_model_4(ch, t) = AIC_values_temp{idx}(4);

    BIC_values_model_1(ch, t) = BIC_values_temp{idx}(1);
    BIC_values_model_2(ch, t) = BIC_values_temp{idx}(2);
    BIC_values_model_3(ch, t) = BIC_values_temp{idx}(3);
    BIC_values_model_4(ch, t) = BIC_values_temp{idx}(4);
end

%% Calculate mean AIC and BIC for each model across all time-point and channel combinations
mean_AIC_model_1 = mean(AIC_values_model_1(:), 'omitnan');
mean_BIC_model_1 = mean(BIC_values_model_1(:), 'omitnan');

mean_AIC_model_2 = mean(AIC_values_model_2(:), 'omitnan');
mean_BIC_model_2 = mean(BIC_values_model_2(:), 'omitnan');

mean_AIC_model_3 = mean(AIC_values_model_3(:), 'omitnan');
mean_BIC_model_3 = mean(BIC_values_model_3(:), 'omitnan');

mean_AIC_model_4 = mean(AIC_values_model_4(:), 'omitnan');
mean_BIC_model_4 = mean(BIC_values_model_4(:), 'omitnan');

%% mean AIC and BIC results for all models
comparison_table = table({'Model 1: AgeGroup + Condition + MeanRT + Direction'; 'Model 2: AgeGroup * Condition + MeanRT + Direction'; 'Model 3: AgeGroup * Condition * MeanRT + Direction'; 'Model 4: AgeGroup + Condition * MeanRT + Direction'}, ...
    [mean_AIC_model_1; mean_AIC_model_2; mean_AIC_model_3; mean_AIC_model_4], ...
    [mean_BIC_model_1; mean_BIC_model_2; mean_BIC_model_3; mean_BIC_model_4], ...
    'VariableNames', {'Model', 'Mean_AIC', 'Mean_BIC'});

disp(comparison_table);

%% Statistical significance of model comparisons

% Threshold for significance
alpha = 0.05;

% Count significant p-values for each comparison
significant_1_vs_2 = sum(p_values_1_vs_2(:) < alpha, 'omitnan');
significant_2_vs_3 = sum(p_values_2_vs_3(:) < alpha, 'omitnan');
significant_3_vs_4 = sum(p_values_3_vs_4(:) < alpha, 'omitnan');
significant_1_vs_4 = sum(p_values_1_vs_4(:) < alpha, 'omitnan');
significant_1_vs_3 = sum(p_values_1_vs_3(:) < alpha, 'omitnan'); % New comparison
significant_2_vs_4 = sum(p_values_2_vs_4(:) < alpha, 'omitnan'); % New comparison

% Proportion of significant comparisons
total_comparisons = channels * times;
proportion_1_vs_2 = significant_1_vs_2 / total_comparisons;
proportion_2_vs_3 = significant_2_vs_3 / total_comparisons;
proportion_3_vs_4 = significant_3_vs_4 / total_comparisons;
proportion_1_vs_4 = significant_1_vs_4 / total_comparisons;
proportion_1_vs_3 = significant_1_vs_3 / total_comparisons; % New comparison
proportion_2_vs_4 = significant_2_vs_4 / total_comparisons; % New comparison

comparison = {'Model 1 vs 2', 'Model 2 vs 3', 'Model 3 vs 4', 'Model 1 vs 4', 'Model 1 vs 3', 'Model 2 vs 4'};
proportions = [proportion_1_vs_2; proportion_2_vs_3; proportion_3_vs_4; proportion_1_vs_4; proportion_1_vs_3; proportion_2_vs_4];

disp(table(comparison', proportions, 'VariableNames', {'Comparison', 'Proportion_Significant'}));
