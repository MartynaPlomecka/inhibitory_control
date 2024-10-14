% EEG TFCE Analysis with Corrected Permutation for Between-Subject and Within-Subject Variables

restoredefaultpath;
clear
clc

%% setup
functions = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\toolboxes\lmeEEG-main\functions';
addpath(genpath(functions));
dependencies = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\toolboxes\ept_TFCE-matlab-master\TFCE\Dependencies';
addpath(genpath(dependencies));

addpath('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_onbinn_theta\fun')
colormaps_path = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_onbinn_theta\mat_files\';
load([colormaps_path 'colormaps.mat']);
load([colormaps_path 'unicolor_maps.mat']);
create_colormap;
mycolormap = customcolormap_preset('red-white-blue');

eeglabPath = '\\psyger-stor02.d.uzh.ch\methlab\Students\Marius\toolboxes\eeglab2021.1';
if ~exist('pop_reref', 'file')
    addpath(eeglabPath);
    eeglab; close;
end

HBNpath = '\\psyger-stor02.d.uzh.ch\methlab\HBN\EEG\Alpha analysis\FOOOF_HBN_final';
load([HBNpath filesep 'chanlocs105.mat']);
chanlocs = chanlocs105;

%% data
addpath('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_on_erps\saccade_locked')
load('all_subjects_sacclocked_ERP_lme_table_separated_direction.mat')

all_subjects = all_subject_tables.SubjectID;
all_ag_rep = all_subject_tables.AgeGroup;
all_cond = all_subject_tables.Condition;
all_dir = all_subject_tables.Direction;
all_mean_rt = all_subject_tables.MeanRT;
all_eeg_data = all_subject_tables.EEGData;

% Convert data to matrix form
n_samples = length(all_subjects);
[channels, times] = size(all_eeg_data{1});
time_points = linspace(-800, 798, times);
all_eeg = cat(3, all_eeg_data{:});
all_eeg = permute(all_eeg, [3, 1, 2]);

% Data normalization and time window selection
time_window = [-500, 100]; % ms
window_indices = time_points >= time_window(1) & time_points <= time_window(2);
all_eeg = all_eeg(:, :, window_indices);
times = sum(window_indices);
time_points = time_points(window_indices);

% Z-score normalization
all_eeg_normalized = (all_eeg - mean(all_eeg(:))) / std(all_eeg(:));
all_mean_rt_normalized = (all_mean_rt - mean(all_mean_rt)) / std(all_mean_rt);

AgeGroup = categorical(all_ag_rep, [0, 1], {'young', 'old'});
Condition = categorical(all_cond, {'pro', 'anti'}, 'Ordinal', false);
Direction = categorical(all_dir, {'left', 'right'}, 'Ordinal', false);
SubjectID = categorical(all_subjects);

% Prep for modeling
rep_eeg_data = squeeze(all_eeg(:, 1, 1));
rep_eeg_tbl = table(rep_eeg_data, AgeGroup, SubjectID, Condition, Direction, all_mean_rt_normalized, ...
    'VariableNames', {'EEG', 'AgeGroup', 'SubjectID', 'Condition', 'Direction', 'MeanRT'});

% Fit initial linear mixed-effects model to extract the design matrix
lme_model = fitlme(rep_eeg_tbl, 'EEG ~ AgeGroup * Condition + Direction + MeanRT + (1|SubjectID)', ...
    'FitMethod', 'ML', 'DummyVarCoding', 'effects');
X = designMatrix(lme_model);

% Perform mass univariate linear regressions
t_obs = nan(channels, times, size(X, 2));
betas = nan(channels, times, size(X, 2));
se = nan(channels, times, size(X, 2));

disp(' mass univariate regressions')
parfor ch = 1:channels
    for t = 1:times
        eeg = squeeze(all_eeg_normalized(:, ch, t));
        [t_obs(ch, t, :), betas(ch, t, :), se(ch, t, :)] = lmeEEG_regress(eeg, X);
    end
end

%% permutation:
nperms = 1000;
ChN = ept_ChN2(chanlocs);

% Generate permutations considering bin-level MeanRT and shuffling AgeGroup
[permutedData] = lmeEEG_permutations_bin_level(nperms, SubjectID, AgeGroup, Condition, Direction, all_mean_rt_normalized);

maxTFCE = nan(nperms, size(X, 2));
disp('Computing permutations...')

placeholder_EEG = rep_eeg_data;

parfor p = 1:nperms
    t_perms = nan(channels, times, size(X, 2));

    % Get permuted labels
    permutedAgeGroup = permutedData{p}.AgeGroup;
    permutedCondition = permutedData{p}.Condition;
    permutedDirection = permutedData{p}.Direction;
    permutedMeanRT = permutedData{p}.MeanRT;

    % Create permuted design matrix
    permuted_tbl = table(placeholder_EEG, permutedAgeGroup, SubjectID, permutedCondition, permutedDirection, permutedMeanRT, ...
        'VariableNames', {'EEG', 'AgeGroup', 'SubjectID', 'Condition', 'Direction', 'MeanRT'});

    % Recompute the design matrix with permuted labels
    lme_model_perm = fitlme(permuted_tbl, 'EEG ~ AgeGroup * Condition + Direction + MeanRT + (1|SubjectID)', ...
        'FitMethod', 'ML', 'DummyVarCoding', 'effects');
    X_perm = designMatrix(lme_model_perm);

    for ch = 1:channels
        EEG_ch = squeeze(all_eeg_normalized(:, ch, :));
        t_perms_ch = nan(times, size(X_perm, 2));
        for t = 1:times
            EEG = EEG_ch(:, t);
            [t_perms_ch(t, :)] = lmeEEG_regress(EEG, X_perm);
        end
        t_perms(ch, :, :) = t_perms_ch;
    end
    maxTFCE_p = nan(1, size(X_perm, 2));
    for i = 2:size(X_perm, 2)
        TFCE_Perm = ept_mex_TFCE2D(t_perms(:, :, i), ChN, [0.66 2]);
        maxTFCE_p(i) = max(abs(TFCE_Perm(:)));
    end
    maxTFCE(p, :) = maxTFCE_p;
end

disp('Permutations completed');

% TFCE results
Results = struct();
for i = 2:size(X, 2)
    coeffName = matlab.lang.makeValidName(lme_model.CoefficientNames{i});
    Results.(coeffName) = lmeEEG_TFCE(squeeze(t_obs(:, :, i)), maxTFCE(:, i), chanlocs, [0.66 2]);
end

save('sacclock_ERP_lme_TFCE_Results.mat', 'Results', 'betas', 't_obs', 'X');

%% Electrode sorting for visualization
T = struct2table(chanlocs);
[~, ~, Th, Rd, ~] = readlocs(chanlocs);
Th = pi / 180 * Th;
[x_coords, y_coords] = pol2cart(Th, Rd);
T.x = x_coords';
T.y = y_coords';
[sortedT, sort_idx] = sortrows(T, 'x', 'descend');

frontal = sortedT(sortedT.x > 0.24, :);
f_idx = sort_idx(sortedT.x > 0.24);
cen_par = sortedT(sortedT.x <= 0.24 & sortedT.x > -0.31, :);
cp_idx = sort_idx(sortedT.x <= 0.24 & sortedT.x > -0.31);
occ = sortedT(sortedT.x <= -0.31, :);
o_idx = sort_idx(sortedT.x <= -0.31);

[~, i1] = sortrows(frontal, 'y', 'ascend');
[~, i2] = sortrows(cen_par, 'y', 'ascend');
[~, i3] = sortrows(occ, 'y', 'ascend');

f_idx = f_idx(i1);
cp_idx = cp_idx(i2);
o_idx = o_idx(i3);

final_idx = [f_idx; cp_idx; o_idx];

% Rearrange t_obs
t_obs_sorted = t_obs(final_idx, :, :);

%% Visualization
% Visualization 1: Significant Regions and White Non-Significant Values
f1 = figure('Color', 'w', 'Position', [5, 10, 2000, 600]);

numCoefficients = length(lme_model.CoefficientNames);
coefficientsToPlot = 2:numCoefficients;
numPlots = length(coefficientsToPlot);

frontal_length = length(f_idx);
cp_length = length(cp_idx);
occ_length = length(o_idx);

subplotMargin = 0.05;
subplotWidth = (1 - (numPlots + 1) * subplotMargin) / numPlots;

for idx = 1:numPlots
    i = coefficientsToPlot(idx);
    subplot('Position', [(idx-1)*(subplotWidth + subplotMargin) + subplotMargin, 0.2, subplotWidth, 0.6]);

    coeffName = matlab.lang.makeValidName(lme_model.CoefficientNames{i});
    significant_mask = Results.(coeffName).Mask;

    t_obs_coeff_sorted = t_obs_sorted(:, :, i);
    significant_mask_sorted = significant_mask(final_idx, :);

    rgb_image = ind2rgb(gray2ind(mat2gray(t_obs_coeff_sorted, [-6 6]), 256), mycolormap);

    % Set non-significant values to white
    for j = 1:3
        temp = rgb_image(:,:,j);
        temp(~significant_mask_sorted) = 1;  
        rgb_image(:,:,j) = temp;
    end

    image(time_points, 1:channels, rgb_image);
    set(gca, 'YDir', 'normal');

    xlabel('Time [ms]', 'FontSize', 10, 'Interpreter', 'none');

    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);

    title(strrep(lme_model.CoefficientNames{i}, '_', ' '), 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');

    if idx == numPlots
        c = colorbar;
        c.Position = [c.Position(1) + 0.05, c.Position(2), c.Position(3)*0.8, c.Position(4)];
        clim([-6 6]);  
    end

    hold on;
    xline(0, 'k', 'LineWidth', 1);
    hold off;

    xlim(time_window);

    h1 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length/2, 'Frontal', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h2 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length/2, ['Centro-' newline 'parietal'], 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h3 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length + occ_length/2, 'Occipital', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');

    set(h1, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h2, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h3, 'Rotation', 35, 'HorizontalAlignment', 'right');
end

colormap(mycolormap);

%% Visualization of the masked version transparent
f2 = figure('Color', 'w', 'Position', [5, 10, 2000, 600]);

for idx = 1:numPlots
    i = coefficientsToPlot(idx);
    subplot('Position', [(idx-1)*(subplotWidth + subplotMargin) + subplotMargin, 0.2, subplotWidth, 0.6]);

    coeffName = matlab.lang.makeValidName(lme_model.CoefficientNames{i});
    significant_mask = Results.(coeffName).Mask;

    t_obs_coeff_sorted = t_obs_sorted(:, :, i);

    significant_mask_sorted = significant_mask(final_idx, :);

    alpha_mask = 0.25 * ~significant_mask_sorted + 1.0 * significant_mask_sorted;

    imagesc(time_points, 1:channels, t_obs_coeff_sorted, 'AlphaData', alpha_mask);

    set(gca, 'YDir', 'normal');
    clim([-6 6]);
    xlabel('Time [ms]', 'FontSize', 10, 'Interpreter', 'none');

    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);

    title(strrep(lme_model.CoefficientNames{i}, '_', ' '), 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');

    if idx == numPlots
        c = colorbar;
        c.Position = [c.Position(1) + 0.05, c.Position(2), c.Position(3)*0.8, c.Position(4)];
    end

    colormap(mycolormap);

    hold on;
    xline(0, 'k', 'LineWidth', 1);
    hold off;

    xlim(time_window);

    h1 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length/2, 'Frontal', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h2 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length/2, ['Centro-' newline 'parietal'], 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h3 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length + occ_length/2, 'Occipital', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');

    set(h1, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h2, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h3, 'Rotation', 35, 'HorizontalAlignment', 'right');
end

%% Visualization of uncorrected data
f3 = figure('Color', 'w', 'Position', [5, 10, 2000, 600]);

for idx = 1:numPlots
    i = coefficientsToPlot(idx);
    subplot('Position', [(idx-1)*(subplotWidth + subplotMargin) + subplotMargin, 0.2, subplotWidth, 0.6]);

    coeffName = matlab.lang.makeValidName(lme_model.CoefficientNames{i});
    t_obs_coeff_sorted = t_obs_sorted(:, :, i);

    imagesc(time_points, 1:channels, t_obs_coeff_sorted);
    set(gca, 'YDir', 'normal');
    clim([-6 6]);
    xlabel('Time [ms]', 'FontSize', 10, 'Interpreter', 'none');

    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);

    title(strrep(lme_model.CoefficientNames{i}, '_', ' '), 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');

    if idx == numPlots
        c = colorbar;
        c.Position = [c.Position(1) + 0.05, c.Position(2), c.Position(3)*0.8, c.Position(4)];
    end

    colormap(mycolormap);

    hold on;
    xline(0, 'k', 'LineWidth', 1);
    hold off;

    xlim(time_window);

    h1 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length/2, 'Frontal', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h2 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length/2, ['Centro-' newline 'parietal'], 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');
    h3 = text(time_window(1) - (time_window(2)-time_window(1))/10, frontal_length + cp_length + occ_length/2, 'Occipital', 'FontName', 'Segoe UI', 'FontSize', 10, 'FontWeight', 'normal');

    set(h1, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h2, 'Rotation', 35, 'HorizontalAlignment', 'right');
    set(h3, 'Rotation', 35, 'HorizontalAlignment', 'right');
end

%% Extract significant electrodes and time windows for age*condition interaction
interaction_coeff = 'AgeGroup_young_Condition_pro';

interaction_mask = Results.(interaction_coeff).Mask;

original_electrode_labels = {chanlocs.labels};

%  mapping between electrode labels and their numbers
electrode_numbers = 1:length(chanlocs);
electrode_map = containers.Map({chanlocs.labels}, num2cell(electrode_numbers));

[sig_electrodes, sig_times] = find(interaction_mask);

sig_electrodes_original = original_electrode_labels(sig_electrodes);

unique_electrodes = unique(sig_electrodes_original);
electrode_time_ranges = cell(length(unique_electrodes), 1);


for i = 1:length(unique_electrodes)
    electrode = unique_electrodes{i};
    electrode_times = time_points(sig_times(sig_electrodes == find(strcmp(original_electrode_labels, electrode))));
    time_ranges = findContinuousacceRanges(electrode_times);
    electrode_time_ranges{i} = time_ranges;

 
end

disp('Significant electrodes and time windows for age*condition interaction:');
disp(['Interaction term: ' interaction_coeff]);
disp('All time points are relative to saccade onset (0 ms)');

for i = 1:length(unique_electrodes)
    electrode = unique_electrodes{i};
    time_ranges = electrode_time_ranges{i};
    electrode_number = electrode_map(electrode);

    disp(['Electrode ', electrode, ' (', num2str(electrode_number), '):']);
    for j = 1:size(time_ranges, 1)
        disp(['  ', num2str(time_ranges(j, 1)), ' to ', num2str(time_ranges(j, 2)), ' ms']);
    end
end

%%
disp('Significant periods between -300 and 0 ms after sacculus onset:');
for i = 1:length(unique_electrodes)
    electrode = unique_electrodes{i};
    time_ranges = electrode_time_ranges{i};
    electrode_number = electrode_map(electrode);

    valid_ranges = time_ranges(time_ranges(:,1) >= -300 & time_ranges(:,2) <=0, :);

    if ~isempty(valid_ranges)
        disp(['Electrode ', electrode, ' (', num2str(electrode_number), '):']);
        for j = 1:size(valid_ranges, 1)
            disp(['  ', num2str(valid_ranges(j, 1)), ' to ', num2str(valid_ranges(j, 2)), ' ms']);
        end
    end
end

frontal_electrodes = sum(ismember(unique_electrodes, frontal.labels));
cp_electrodes = sum(ismember(unique_electrodes, cen_par.labels));
occ_electrodes = sum(ismember(unique_electrodes, occ.labels));

disp('Number of significant electrodes by brain region:');
disp(['  Frontal: ' num2str(frontal_electrodes)]);
disp(['  Centro-parietal: ' num2str(cp_electrodes)]);
disp(['  Occipital: ' num2str(occ_electrodes)]);


%% Functions

function [tval, b, se] = lmeEEG_regress(y, X)
[n, ncolX] = size(X);
[Q, R, perm] = qr(X, 0);
if isempty(R)
    p = 0;
elseif isvector(R)
    p = double(abs(R(1)) > 0);
else
    p = sum(abs(diag(R)) > max(n, ncolX) * eps(R(1)));
end
if p < ncolX
    warning('Rank deficient design matrix.');
    R = R(1:p, 1:p);
    Q = Q(:, 1:p);
    perm = perm(1:p);
end
b = zeros(ncolX, 1);
b(perm) = R \ (Q' * y);
RI = R \ eye(p);
nu = max(0, n - p);
yhat = X * b;
r = y - yhat;
normr = norm(r);
if nu ~= 0
    rmse = normr / sqrt(nu);
else
    rmse = NaN;
end
se = zeros(ncolX, 1);
se(perm, :) = rmse * sqrt(sum(abs(RI).^2, 2));
tval = b ./ se;
end

function [Results] = lmeEEG_TFCE(T_Obs, maxTFCE, e_loc, E_H)
if ~isequal(size(T_Obs, 1), length(e_loc))
    error('Number of channels in data does not equal that of locations file')
end
ChN = ept_ChN2(e_loc);
if ismatrix(T_Obs)
    TFCE_Obs = ept_mex_TFCE2D(T_Obs, ChN, E_H);
elseif ndims(T_Obs) == 3
    TFCE_Obs = ept_mex_TFCE3D(T_Obs, ChN, E_H);
end
Alpha = .05;
nPerm = length(maxTFCE);
maxTFCE = sort([maxTFCE; max(abs(TFCE_Obs(:)))]);
maxTFCEcrit = maxTFCE(round(nPerm * (1 - Alpha)));
Mask = abs(TFCE_Obs) >= maxTFCEcrit;
P_Values = NaN(size(TFCE_Obs, 1), size(TFCE_Obs, 2));
for idx = 1:size(TFCE_Obs, 1)
    for jdx = 1:size(TFCE_Obs, 2)
        P_Values(idx, jdx) = sum(abs(TFCE_Obs(idx, jdx)) <= maxTFCE) / (nPerm + 1);
    end
end
Results.Obs = T_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = maxTFCE;
Results.P_Values = P_Values;
Results.Mask = Mask;
end

function [permutedData] = lmeEEG_permutations_bin_level(nperms, SubjectID, AgeGroup, Condition, Direction, MeanRT)
uniqueSubjects = unique(SubjectID);
numSubjects = length(uniqueSubjects);
permutedData = cell(nperms, 1);

for p = 1:nperms
    permutedAgeGroup = categorical(repmat(missing, size(SubjectID)), categories(AgeGroup));
    permutedCondition = categorical(repmat(missing, size(Condition)), categories(Condition));
    permutedDirection = categorical(repmat(missing, size(Direction)), categories(Direction));
    permutedMeanRT = zeros(size(MeanRT));

    shuffledGroups = AgeGroup(randperm(numSubjects));

    for idx = 1:numSubjects
        subj = uniqueSubjects(idx);
        subjIdx = SubjectID == subj;

        permutedAgeGroup(subjIdx) = shuffledGroups(idx);

        bins = table(Condition(subjIdx), Direction(subjIdx), MeanRT(subjIdx), 'VariableNames', {'Condition', 'Direction', 'MeanRT'});

        permIdx = randperm(height(bins));
        permutedCondition(subjIdx) = bins.Condition(permIdx);
        permutedDirection(subjIdx) = bins.Direction(permIdx);
        permutedMeanRT(subjIdx) = bins.MeanRT(permIdx);
    end

    permutedData{p} = struct('AgeGroup', permutedAgeGroup, 'Condition', permutedCondition, 'Direction', permutedDirection, 'MeanRT', permutedMeanRT);
end
end

function ranges = findContinuousacceRanges(times)
if isempty(times)
    ranges = [];
    return;
end

times = sort(times);
diffs = diff(times);
breaks = find(diffs > 2);  % TBC We allow for 2ms gap to be considered continuous

if isempty(breaks)
    ranges = [times(1), times(end)];
else
    ranges = zeros(length(breaks) + 1, 2);
    ranges(1, :) = [times(1), times(breaks(1))];
    for i = 2:length(breaks)
        ranges(i, :) = [times(breaks(i-1) + 1), times(breaks(i))];
    end
    ranges(end, :) = [times(breaks(end) + 1), times(end)];
end
end
