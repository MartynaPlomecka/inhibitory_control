% Define paths for EEGLAB and data
eeglabPath = '\\psyger-stor02.d.uzh.ch\methlab\Students\Marius\toolboxes\eeglab2021.1';
dataPath = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\dataT1_bas\test_baselined_clean_longsegmented_data_correcttrials_only_final_unfoldclean_opticat\';

if ~exist('pop_reref', 'file')
    addpath(eeglabPath);
    eeglab; close;
end

% Define clusters
midfrontal_cluster = [15, 9, 3];
centroparietal_cluster = [45, 46, 66, 51, 52, 65];
occipital_cluster = [58, 62, 69, 61, 68];
clusters = {midfrontal_cluster, centroparietal_cluster, occipital_cluster};
cluster_names = {'Midfrontal', 'Centroparietal', 'Occipital'};

% Create parallel pool if not existing
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 16);
end

% Load data and process per subject
d = dir(dataPath);
d = d([d.isdir]);
d(startsWith({d.name}, '.')) = [];
results = cell(numel(d), 1);

parfor sj = 1:numel(d)
    try
        disp(['Processing subject: ' num2str(sj)]);
        subjectPath = fullfile(d(sj).folder, d(sj).name);
        eegFile = dir(fullfile(subjectPath, '*sacclockedEEG.mat'));
        eegData = load(fullfile(eegFile.folder, eegFile.name));
        saccEEG = eegData.saccEEG;

        % Extract sacculus events and age group
        sacculus_events = saccEEG.event(2:2:end);
        age_group = sacculus_events(1).age; % 0 is young, 1 is old

        % Extract pro and anti trials and RTs
        pro_indices = find(strcmp({sacculus_events.cond}, 'pro'));
        anti_indices = find(strcmp({sacculus_events.cond}, 'anti'));
        pro_rt = [sacculus_events(pro_indices).rt];
        anti_rt = [sacculus_events(anti_indices).rt];

        % Bin trials and get mean RTs for each bin
        [pro_binned, ~, ~, pro_mean_rt] = bin_trials(pro_indices, pro_rt);
        [anti_binned, ~, ~, anti_mean_rt] = bin_trials(anti_indices, anti_rt);

        % Initialize ERPs for each cluster and bin
        pro_erp = cell(3, 3); % {cluster}{bin}
        anti_erp = cell(3, 3);
        for c = 1:3
            for bin = 1:3
                % Pro ERPs
                if ~isempty(pro_binned{bin})
                    pro_bin_indices = ismember([sacculus_events.urevent], [sacculus_events(pro_binned{bin}).urevent]);
                    pro_erp{c, bin} = mean(saccEEG.data(clusters{c}, :, pro_bin_indices), 3);
                else
                    pro_erp{c, bin} = [];
                end
                % Anti ERPs
                if ~isempty(anti_binned{bin})
                    anti_bin_indices = ismember([sacculus_events.urevent], [sacculus_events(anti_binned{bin}).urevent]);
                    anti_erp{c, bin} = mean(saccEEG.data(clusters{c}, :, anti_bin_indices), 3);
                else
                    anti_erp{c, bin} = [];
                end
            end
        end

        results{sj} = struct('age_group', age_group, 'pro_erp', {pro_erp}, 'anti_erp', {anti_erp}, ...
            'times', saccEEG.times, 'pro_mean_rt', pro_mean_rt, 'anti_mean_rt', anti_mean_rt);

    catch ME
        warning(['Error with subject ' num2str(sj) ': ' ME.message]);
        results{sj} = struct('age_group', NaN, 'pro_erp', {cell(3,3)}, 'anti_erp', {cell(3,3)}, 'times', [], ...
            'pro_mean_rt', [], 'anti_mean_rt', []);
    end
end

% Initialize variables to collect ERPs
young_pro_erps = cell(3, 3);
young_anti_erps = cell(3, 3);
old_pro_erps = cell(3, 3);
old_anti_erps = cell(3, 3);
time = [];

% Collect ERPs into young and old groups
for sj = 1:numel(results)
    if ~isnan(results{sj}.age_group)
        if results{sj}.age_group == 0 % Young group
            for c = 1:3
                for bin = 1:3
                    if ~isempty(results{sj}.pro_erp{c, bin})
                        young_pro_erps{c, bin} = cat(3, young_pro_erps{c, bin}, results{sj}.pro_erp{c, bin});
                    end
                    if ~isempty(results{sj}.anti_erp{c, bin})
                        young_anti_erps{c, bin} = cat(3, young_anti_erps{c, bin}, results{sj}.anti_erp{c, bin});
                    end
                end
            end
        elseif results{sj}.age_group == 1 % Old group
            for c = 1:3
                for bin = 1:3
                    if ~isempty(results{sj}.pro_erp{c, bin})
                        old_pro_erps{c, bin} = cat(3, old_pro_erps{c, bin}, results{sj}.pro_erp{c, bin});
                    end
                    if ~isempty(results{sj}.anti_erp{c, bin})
                        old_anti_erps{c, bin} = cat(3, old_anti_erps{c, bin}, results{sj}.anti_erp{c, bin});
                    end
                end
            end
        end
        if isempty(time) && ~isempty(results{sj}.times)
            time = results{sj}.times;
        end
    end
end
% Average and SEM calculations for ERPs
young_pro_avg = cell(3, 3); young_anti_avg = cell(3, 3);
old_pro_avg = cell(3, 3); old_anti_avg = cell(3, 3);
young_pro_sem = cell(3, 3); young_anti_sem = cell(3, 3);
old_pro_sem = cell(3, 3); old_anti_sem = cell(3, 3);

for c = 1:3
    for bin = 1:3
        if ~isempty(young_pro_erps{c, bin})
            young_pro_avg{c, bin} = mean(young_pro_erps{c, bin}, 3);
            young_pro_sem{c, bin} = std(young_pro_erps{c, bin}, 0, 3) / sqrt(size(young_pro_erps{c, bin}, 3));
        end
        if ~isempty(young_anti_erps{c, bin})
            young_anti_avg{c, bin} = mean(young_anti_erps{c, bin}, 3);
            young_anti_sem{c, bin} = std(young_anti_erps{c, bin}, 0, 3) / sqrt(size(young_anti_erps{c, bin}, 3));
        end
        if ~isempty(old_pro_erps{c, bin})
            old_pro_avg{c, bin} = mean(old_pro_erps{c, bin}, 3);
            old_pro_sem{c, bin} = std(old_pro_erps{c, bin}, 0, 3) / sqrt(size(old_pro_erps{c, bin}, 3));
        end
        if ~isempty(old_anti_erps{c, bin})
            old_anti_avg{c, bin} = mean(old_anti_erps{c, bin}, 3);
            old_anti_sem{c, bin} = std(old_anti_erps{c, bin}, 0, 3) / sqrt(size(old_anti_erps{c, bin}, 3));
        end
    end
end

% Calculate mean RT values for plotting - separate by age and condition
young_pro_rt = []; young_anti_rt = [];
old_pro_rt = []; old_anti_rt = [];

for sj = 1:numel(results)
    if isstruct(results{sj}) && ~isnan(results{sj}.age_group)
        if results{sj}.age_group == 0  % Young
            if ~isempty(results{sj}.pro_mean_rt)
                young_pro_rt = [young_pro_rt; results{sj}.pro_mean_rt];
            end
            if ~isempty(results{sj}.anti_mean_rt)
                young_anti_rt = [young_anti_rt; results{sj}.anti_mean_rt];
            end
        else  % Old
            if ~isempty(results{sj}.pro_mean_rt)
                old_pro_rt = [old_pro_rt; results{sj}.pro_mean_rt];
            end
            if ~isempty(results{sj}.anti_mean_rt)
                old_anti_rt = [old_anti_rt; results{sj}.anti_mean_rt];
            end
        end
    end
end

%% Plot preparations
% Calculate means for each group (negative for saccade-locked)
avg_young_pro_rt = -mean(young_pro_rt, 1);
avg_young_anti_rt = -mean(young_anti_rt, 1);
avg_old_pro_rt = -mean(old_pro_rt, 1);
avg_old_anti_rt = -mean(old_anti_rt, 1);

% y-axis limits for each cluster
cluster_ylims = {[-1.29 0.77],    % Midfrontal
    [-0.6 3],         % Centroparietal
    [-1.2 1.2]};      % Occipital

%% Original 2x2 plot
for c = 1:3
    figure('Color', 'white', 'Position', [100, 100, 1200, 800]);
    colors = {[0, 0, 0], [0, 0, 1], [0.7, 0.9, 1]}; % black, blue, light blue
    bin_names = {'Fast', 'Medium', 'Slow'};

    for i = 1:4
        subplot(2, 2, i);
        hold on;
        h = zeros(1, 3);
        ylims = cluster_ylims{c};
        ylim(ylims);

        line([0 0], ylims, 'Color', 'r', 'LineStyle', '-', 'HandleVisibility', 'off');

        for bin = 1:3
            switch i
                case 1 % Young Pro
                    erp_data = young_pro_avg{c, bin};
                    sem_data = young_pro_sem{c, bin};
                    mean_rt = avg_young_pro_rt(bin);
                    title_str = 'Young Pro';
                case 2 % Young Anti
                    erp_data = young_anti_avg{c, bin};
                    sem_data = young_anti_sem{c, bin};
                    mean_rt = avg_young_anti_rt(bin);
                    title_str = 'Young Anti';
                case 3 % Old Pro
                    erp_data = old_pro_avg{c, bin};
                    sem_data = old_pro_sem{c, bin};
                    mean_rt = avg_old_pro_rt(bin);
                    title_str = 'Old Pro';
                case 4 % Old Anti
                    erp_data = old_anti_avg{c, bin};
                    sem_data = old_anti_sem{c, bin};
                    mean_rt = avg_old_anti_rt(bin);
                    title_str = 'Old Anti';
            end

            if ~isempty(erp_data)
                temp = shadedErrorBar(time, mean(erp_data, 1), mean(sem_data, 1), ...
                    {'Color', colors{bin}, 'LineStyle', '-', 'LineWidth', 2.5}, 1);
                h(bin) = temp.mainLine;
                temp.patch.FaceColor = colors{bin};
                temp.patch.FaceAlpha = 0.3;
                temp.patch.EdgeColor = 'none';

                line([mean_rt mean_rt], ylims, 'Color', colors{bin}, ...
                    'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end

        xlabel('Time (ms)', 'FontSize', 14, 'FontName', 'Helvetica');
        ylabel('Amplitude (µV)', 'FontSize', 14, 'FontName', 'Helvetica');
        title([title_str ' '], 'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'normal');
        legend(h, bin_names, 'FontSize', 12, 'FontName', 'Helvetica');
        xlim([-500, 100]);
        set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
        box on;
    end

    sgtitle(['ERPs sacclocked,  ' cluster_names{c} ' Cluster'], 'FontWeight', 'normal');
    saveas(gcf, fullfile(['ERPs_RT_bins_' cluster_names{c} '.png']));
end

%% Pro vs Anti Comparison (All Ages Combined)
% Calculate age-averaged ERPs
pro_avg = cell(3, 3); anti_avg = cell(3, 3);
pro_sem = cell(3, 3); anti_sem = cell(3, 3);

for c = 1:3
    for bin = 1:3
        combined_pro = cat(3, young_pro_erps{c, bin}, old_pro_erps{c, bin});
        if ~isempty(combined_pro)
            pro_avg{c, bin} = mean(combined_pro, 3);
            pro_sem{c, bin} = std(combined_pro, 0, 3) / sqrt(size(combined_pro, 3));
        end

        combined_anti = cat(3, young_anti_erps{c, bin}, old_anti_erps{c, bin});
        if ~isempty(combined_anti)
            anti_avg{c, bin} = mean(combined_anti, 3);
            anti_sem{c, bin} = std(combined_anti, 0, 3) / sqrt(size(combined_anti, 3));
        end
    end
end

% Calculate combined RTs
avg_pro_rt = -mean([young_pro_rt; old_pro_rt], 1);  % Add negative sign
avg_anti_rt = -mean([young_anti_rt; old_anti_rt], 1); % Add negative sign
% Plot Pro vs Anti comparison
for c = 1:3
    figure('Color', 'white', 'Position', [100, 100, 600, 800]);
    pro_color = [0, 0, 1];    % Blue
    anti_color = [1, 0, 0];   % Red

    for bin = 1:3
        subplot(3, 1, bin);
        hold on;
        ylims = cluster_ylims{c};
        ylim(ylims);

        line([0 0], ylims, 'Color', 'k', 'LineStyle', '-', 'HandleVisibility', 'off');

        if ~isempty(pro_avg{c, bin})
            temp_pro = shadedErrorBar(time, mean(pro_avg{c, bin}, 1), mean(pro_sem{c, bin}, 1), ...
                {'Color', pro_color, 'LineStyle', '-', 'LineWidth', 2}, 1);
            temp_pro.patch.FaceColor = pro_color;
            temp_pro.patch.FaceAlpha = 0.2;
            temp_pro.patch.EdgeColor = 'none';

            line([avg_pro_rt(bin) avg_pro_rt(bin)], ylims, 'Color', pro_color, ...
                'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

        if ~isempty(anti_avg{c, bin})
            temp_anti = shadedErrorBar(time, mean(anti_avg{c, bin}, 1), mean(anti_sem{c, bin}, 1), ...
                {'Color', anti_color, 'LineStyle', '-', 'LineWidth', 2}, 1);
            temp_anti.patch.FaceColor = anti_color;
            temp_anti.patch.FaceAlpha = 0.2;
            temp_anti.patch.EdgeColor = 'none';

            line([avg_anti_rt(bin) avg_anti_rt(bin)], ylims, 'Color', anti_color, ...
                'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

        xlabel('Time (ms)', 'FontSize', 12, 'FontName', 'Helvetica');
        ylabel('Amplitude (µV)', 'FontSize', 12, 'FontName', 'Helvetica');

        switch bin
            case 1; title_str = 'Fast Trials';
            case 2; title_str = 'Medium Trials';
            case 3; title_str = 'Slow Trials';
        end
        title(title_str, 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'normal');

        if bin == 1
            legend([temp_pro.mainLine, temp_anti.mainLine], {'Pro', 'Anti'}, ...
                'FontSize', 12, 'FontName', 'Helvetica', 'Location', 'best');
        end

        xlim([-500, 100]);
        set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
        box on;
    end

    sgtitle([cluster_names{c} ' Cluster'], 'FontWeight', 'normal', 'FontSize', 16);
    saveas(gcf, fullfile(['ERPs_RT_bins_ProAnti_' cluster_names{c} '_SaccadeLocked.png']));
end

%% Age Comparison (Pro & Anti Combined)
for c = 1:3
    figure('Color', 'white', 'Position', [100, 100, 600, 800]);
    young_color = [0, 0, 1];    % Blue
    old_color = [1, 0, 0];      % Red

    for bin = 1:3
        subplot(3, 1, bin);
        hold on;
        ylims = cluster_ylims{c};
        ylim(ylims);

        line([0 0], ylims, 'Color', 'k', 'LineStyle', '-', 'HandleVisibility', 'off');

        if ~isempty(young_pro_avg{c, bin}) && ~isempty(young_anti_avg{c, bin})
            young_avg = mean(cat(3, young_pro_avg{c, bin}, young_anti_avg{c, bin}), 3);
            young_sem = std(cat(3, young_pro_avg{c, bin}, young_anti_avg{c, bin}), 0, 3) / ...
                sqrt(size(cat(3, young_pro_avg{c, bin}, young_anti_avg{c, bin}), 3));

            temp_young = shadedErrorBar(time, mean(young_avg, 1), mean(young_sem, 1), ...
                {'Color', young_color, 'LineStyle', '-', 'LineWidth', 2}, 1);
            temp_young.patch.FaceColor = young_color;
            temp_young.patch.FaceAlpha = 0.2;
            temp_young.patch.EdgeColor = 'none';

            young_rt = mean([avg_young_pro_rt(bin), avg_young_anti_rt(bin)]);  % Already negative
            line([young_rt young_rt], ylims, 'Color', young_color, ...
                'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

        if ~isempty(old_pro_avg{c, bin}) && ~isempty(old_anti_avg{c, bin})
            old_avg = mean(cat(3, old_pro_avg{c, bin}, old_anti_avg{c, bin}), 3);
            old_sem = std(cat(3, old_pro_avg{c, bin}, old_anti_avg{c, bin}), 0, 3) / ...
                sqrt(size(cat(3, old_pro_avg{c, bin}, old_anti_avg{c, bin}), 3));

            temp_old = shadedErrorBar(time, mean(old_avg, 1), mean(old_sem, 1), ...
                {'Color', old_color, 'LineStyle', '-', 'LineWidth', 2}, 1);
            temp_old.patch.FaceColor = old_color;
            temp_old.patch.FaceAlpha = 0.2;
            temp_old.patch.EdgeColor = 'none';

            old_rt = mean([avg_old_pro_rt(bin), avg_old_anti_rt(bin)]);  % Already negative
            line([old_rt old_rt], ylims, 'Color', old_color, ...
                'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

        xlabel('Time (ms)', 'FontSize', 12, 'FontName', 'Helvetica');
        ylabel('Amplitude (µV)', 'FontSize', 12, 'FontName', 'Helvetica');

        switch bin
            case 1; title_str = 'Fast Trials';
            case 2; title_str = 'Medium Trials';
            case 3; title_str = 'Slow Trials';
        end
        title(title_str, 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'normal');

        if bin == 1
            legend([temp_young.mainLine, temp_old.mainLine], {'Young', 'Old'}, ...
                'FontSize', 12, 'FontName', 'Helvetica', 'Location', 'best');
        end

        xlim([-500, 100]);
        set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
        box on;
    end

    sgtitle([cluster_names{c} ' Cluster'], 'FontWeight', 'normal', 'FontSize', 16);
    saveas(gcf, fullfile(['ERPs_RT_bins_AgeComparison_' cluster_names{c} '_SaccadeLocked.png']));
end

function [binned_indices, bin_order, sorted_rt, mean_rt] = bin_trials(indices, rt)
[sorted_rt, sort_idx] = sort(rt, 'descend');  % descending order
sorted_indices = indices(sort_idx);
n_trials = length(sorted_indices);
bin_size = floor(n_trials / 3);

% Reverse the bin assignment to match slow->medium->fast
binned_indices = {sorted_indices(2*bin_size+1:end), ...  % fast trials
    sorted_indices(bin_size+1:2*bin_size), ... % medium trials
    sorted_indices(1:bin_size)};  % slow trials

mean_rt = [mean(rt(sort_idx(2*bin_size+1:end))), ...  % fast mean
    mean(rt(sort_idx(bin_size+1:2*bin_size))), ... % medium mean
    mean(rt(sort_idx(1:bin_size)))];  % slow mean

bin_order = {'fast', 'medium', 'slow'};
end