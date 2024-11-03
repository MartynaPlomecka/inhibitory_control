clear all; close all; clc;

eeglabPath = '\\psyger-stor02.d.uzh.ch\methlab\Students\Marius\toolboxes\eeglab2021.1';
dataPath = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\dataT1_bas\test_baselined_clean_longsegmented_data_correcttrials_only_final_unfoldclean_opticat\';
addpath('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_onbinn_theta\fun')
colormaps_path = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_onbinn_theta\mat_files\';
load([colormaps_path 'colormaps.mat']);
load([colormaps_path 'unicolor_maps.mat']);
create_colormap;
mycolormap = customcolormap_preset('red-white-blue');

if ~exist('pop_reref', 'file')
    addpath(eeglabPath);
    eeglab; close;
end

midfrontal_cluster = [15, 9, 3];
centroparietal_cluster = [45, 46, 66, 51, 52, 65];
occipital_cluster = [58, 62, 69, 61, 68];

clusters = {midfrontal_cluster, centroparietal_cluster, occipital_cluster};
cluster_names = {'Midfrontal', 'Centroparietal', 'Occipital'};

% time points for each cluster (in ms)
cluster_times = {
    [150, 250, 350],    % Midfrontal
    [200, 300, 350],    % Centroparietal
    [150, 250, 350]     % Occipital
};

% Time window for averaging (plusminus 10ms)
time_window = 10;

d = dir(dataPath);
d = d([d.isdir]);
d(startsWith({d.name}, '.')) = [];

pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 16);
end

results = cell(numel(d), 1);

parfor sj = 1:numel(d)
    try
        disp(['Processing subject: ' num2str(sj)]);
        
        % Load subject data
        subjectPath = fullfile(d(sj).folder, d(sj).name);
        eegFile = dir(fullfile(subjectPath, '*stimlockedEEG.mat'));
        eegData = load(fullfile(eegFile.folder, eegFile.name));
        stimEEG = eegData.stimEEG;
        
        % Get stimulus events and age group
        stimulus_events = stimEEG.event(1:2:end);
        age_group = stimulus_events(1).age; % 0 is young, 1 is old
        
        % Find pro and anti trials
        pro_indices = find(strcmp({stimulus_events.cond}, 'pro'));
        anti_indices = find(strcmp({stimulus_events.cond}, 'anti'));
        
        % Calculate ERPs for all channels
        pro_erp = mean(stimEEG.data(:, :, pro_indices), 3);
        anti_erp = mean(stimEEG.data(:, :, anti_indices), 3);
        
        results{sj} = struct('age_group', age_group, 'pro_erp', pro_erp, ...
            'anti_erp', anti_erp, 'times', stimEEG.times, 'chanlocs', stimEEG.chanlocs);
        
    catch ME
        warning(['Error with subject ' num2str(sj) ': ' ME.message]);
        results{sj} = struct('age_group', NaN, 'pro_erp', [], ...
            'anti_erp', [], 'times', [], 'chanlocs', []);
    end
end

young_pro_erps = [];
young_anti_erps = [];
old_pro_erps = [];
old_anti_erps = [];
time = [];
chanlocs = [];

n_young = 0;
n_old = 0;

for sj = 1:numel(results)
    if ~isnan(results{sj}.age_group)
        if results{sj}.age_group == 0 % Young group
            if isempty(young_pro_erps)
                young_pro_erps = results{sj}.pro_erp;
                young_anti_erps = results{sj}.anti_erp;
                n_young = 1;
            else
                young_pro_erps = young_pro_erps + results{sj}.pro_erp;
                young_anti_erps = young_anti_erps + results{sj}.anti_erp;
                n_young = n_young + 1;
            end
        elseif results{sj}.age_group == 1 % Old group
            if isempty(old_pro_erps)
                old_pro_erps = results{sj}.pro_erp;
                old_anti_erps = results{sj}.anti_erp;
                n_old = 1;
            else
                old_pro_erps = old_pro_erps + results{sj}.pro_erp;
                old_anti_erps = old_anti_erps + results{sj}.anti_erp;
                n_old = n_old + 1;
            end
        end
        if isempty(time) && ~isempty(results{sj}.times)
            time = results{sj}.times;
            chanlocs = results{sj}.chanlocs;
        end
    end
end

% Calculate averages
young_pro_erps = young_pro_erps / n_young;
young_anti_erps = young_anti_erps / n_young;
old_pro_erps = old_pro_erps / n_old;
old_anti_erps = old_anti_erps / n_old;

% Function to find closest time index
find_time_idx = @(target_time, times) find(abs(times - target_time) == min(abs(times - target_time)));

mycolormap = customcolormap_preset('red-white-blue');

%% Plot

for c = 1:3
    % Create figure with 4 rows (conditions) x 3 columns (time points)
    figure('Color', 'w', 'InvertHardcopy', 'off', 'Position', [100 100 1200 1000]);
    
    all_data = [];
    for t = 1:length(cluster_times{c})
        target_time = cluster_times{c}(t);
        window_start = find_time_idx(target_time - time_window, time);
        window_end = find_time_idx(target_time + time_window, time);
        
        all_data = cat(2, all_data, ...
            mean(young_pro_erps(:, window_start:window_end), 2), ...
            mean(young_anti_erps(:, window_start:window_end), 2), ...
            mean(old_pro_erps(:, window_start:window_end), 2), ...
            mean(old_anti_erps(:, window_start:window_end), 2));
    end
    clim_val = max(abs([min(all_data(:)) max(all_data(:))]));
    
    sgtitle([cluster_names{c} ' Cluster Analysis'], ...
        'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica');  
    
    % Condition labels to the left of each row
    conditions = {'Young Pro', 'Young Anti', 'Old Pro', 'Old Anti'};
    
% Update timing labels below each column
for t = 1:length(cluster_times{c})
    target_time = cluster_times{c}(t);
    annotation('textbox', [0.12 + (t-1)*0.27, 0.01, 0.2, 0.04], ...
        'String', [num2str(target_time) ' ms ± 10 ms'], ...
        'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none');
end

    % Plot each row and condition
    for row = 1:4
        for t = 1:length(cluster_times{c})
            target_time = cluster_times{c}(t);
            window_start = find_time_idx(target_time - time_window, time);
            window_end = find_time_idx(target_time + time_window, time);
            
            % Averaged brain activities for each condition
            data_conditions = {young_pro_erps, young_anti_erps, old_pro_erps, old_anti_erps};
            data = mean(data_conditions{row}(:, window_start:window_end), 2);
            
            % Plot for each time point column
            ax = subplot(4, 3, (row - 1) * 3 + t);
            topoplot(data, chanlocs, 'electrodes', 'on', 'style', 'both', 'conv', 'on');
            colormap(mycolormap);
            colorbar;
            clim([-clim_val clim_val]);
            
            % Add condition label only in the first column of each row
            if t == 1
                text(ax, -0.3, 0.5, conditions{row}, 'FontSize', 16, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                    'Units', 'normalized', 'Rotation', 0);
            end
        end
    end
    
    saveas(gcf, fullfile(['Topoplots_' cluster_names{c} '_cluster.png']));
end
