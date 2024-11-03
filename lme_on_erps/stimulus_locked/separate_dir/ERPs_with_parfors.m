% ERP Analysis Script 
cd \\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\martyna\scripts_for_running_tfr_analysis\lme_on_erps\stimulus_locked\separate_dir
eeglabPath = '\\psyger-stor02.d.uzh.ch\methlab\Students\Marius\toolboxes\eeglab2021.1';
dataPath = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\dataT1_bas\test_baselined_clean_longsegmented_data_correcttrials_only_final_unfoldclean_opticat\';

if ~exist('pop_reref', 'file')
    addpath(eeglabPath);
    eeglab; close;
end

%%
d = dir(dataPath);
d = d([d.isdir]);
d(startsWith({d.name}, '.')) = [];

%E19, E11, E4
midfrontal_cluster = [15, 9, 3];

% 54, 55, 79, 61, 62, 78
centroparietal_cluster = [45, 46, 66, 51, 52, 65];

%70, 75, 83, 74, 82
occipital_cluster = [58, 62, 69, 61, 68];

clusters = {midfrontal_cluster, centroparietal_cluster, occipital_cluster};
cluster_names = {'Midfrontal', 'Centroparietal', 'Occipital'};

pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 16);  
end

results = cell(numel(d), 1);

parfor sj = 1:numel(d)
    try
        disp(['Processing subject: ' num2str(sj)]);
        
        subjectPath = fullfile(d(sj).folder, d(sj).name);
        eegFile = dir(fullfile(subjectPath, '*stimlockedEEG.mat'));
        eegData = load(fullfile(eegFile.folder, eegFile.name));
        stimEEG = eegData.stimEEG;
        
        %stimulus events and age group
        stimulus_events = stimEEG.event(1:2:end);
        age_group = stimulus_events(1).age; % 0 is young, 1 is old
        
        %  pro and anti trials
        pro_indices = find(strcmp({stimulus_events.cond}, 'pro'));
        anti_indices = find(strcmp({stimulus_events.cond}, 'anti'));
        
        % ERPs for this subject (for each cluster)
        pro_erp = cell(1, 3);
        anti_erp = cell(1, 3);
        for c = 1:3
            pro_erp{c} = mean(stimEEG.data(clusters{c}, :, pro_indices), 3);
            anti_erp{c} = mean(stimEEG.data(clusters{c}, :, anti_indices), 3);
        end
        
        results{sj} = struct('age_group', age_group, 'pro_erp', {pro_erp}, 'anti_erp', {anti_erp}, 'times', stimEEG.times);
        
    catch ME
        warning(['Error with subject ' num2str(sj) ': ' ME.message]);
        results{sj} = struct('age_group', NaN, 'pro_erp', {cell(1,3)}, 'anti_erp', {cell(1,3)}, 'times', []);
    end
end

young_pro_erps = cell(1, 3);
young_anti_erps = cell(1, 3);
old_pro_erps = cell(1, 3);
old_anti_erps = cell(1, 3);
time = [];

for c = 1:3
    young_pro_erps{c} = [];
    young_anti_erps{c} = [];
    old_pro_erps{c} = [];
    old_anti_erps{c} = [];
end

for sj = 1:numel(results)
    if ~isnan(results{sj}.age_group)
        if results{sj}.age_group == 0 % Young group
            for c = 1:3
                young_pro_erps{c} = cat(3, young_pro_erps{c}, results{sj}.pro_erp{c});
                young_anti_erps{c} = cat(3, young_anti_erps{c}, results{sj}.anti_erp{c});
            end
        elseif results{sj}.age_group == 1 % Old group
            for c = 1:3
                old_pro_erps{c} = cat(3, old_pro_erps{c}, results{sj}.pro_erp{c});
                old_anti_erps{c} = cat(3, old_anti_erps{c}, results{sj}.anti_erp{c});
            end
        end
        if isempty(time) && ~isempty(results{sj}.times)
            time = results{sj}.times;
        end
    end
end

% Grand average ERPs and SEMs for each cluster
young_pro_avg = cell(1, 3);
young_anti_avg = cell(1, 3);
old_pro_avg = cell(1, 3);
old_anti_avg = cell(1, 3);
young_pro_sem = cell(1, 3);
young_anti_sem = cell(1, 3);
old_pro_sem = cell(1, 3);
old_anti_sem = cell(1, 3);

for c = 1:3
    % GA+ ERPs
    young_pro_avg{c} = mean(young_pro_erps{c}, 3);
    young_anti_avg{c} = mean(young_anti_erps{c}, 3);
    old_pro_avg{c} = mean(old_pro_erps{c}, 3);
    old_anti_avg{c} = mean(old_anti_erps{c}, 3);

    young_pro_sem{c} = std(young_pro_erps{c}, 0, 3) / sqrt(size(young_pro_erps{c}, 3));
    young_anti_sem{c} = std(young_anti_erps{c}, 0, 3) / sqrt(size(young_anti_erps{c}, 3));
    old_pro_sem{c} = std(old_pro_erps{c}, 0, 3) / sqrt(size(old_pro_erps{c}, 3));
    old_anti_sem{c} = std(old_anti_erps{c}, 0, 3) / sqrt(size(old_anti_erps{c}, 3));
end

%%
%plots for each cluster
for c = 1:3
    figure('Color', 'white'); 
    hold on;

    h1 = shadedErrorBar(time, mean(young_pro_avg{c}, 1), mean(young_pro_sem{c}, 1), {'Color', [0.6, 0.6, 1], 'LineStyle', '-', 'LineWidth', 2});
    h2 = shadedErrorBar(time, mean(young_anti_avg{c}, 1), mean(young_anti_sem{c}, 1), {'Color', [1, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 2});
    h3 = shadedErrorBar(time, mean(old_pro_avg{c}, 1), mean(old_pro_sem{c}, 1), {'Color', [0, 0, 1], 'LineStyle', '--', 'LineWidth', 2});
    h4 = shadedErrorBar(time, mean(old_anti_avg{c}, 1), mean(old_anti_sem{c}, 1), {'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 2});

    xlabel('Time (ms)', 'FontSize', 14, 'FontName', 'Helvetica');
    ylabel('Amplitude (µV)', 'FontSize', 14, 'FontName', 'Helvetica');
    title(['STIMULUS LOCKED ERP: Pro vs Anti Saccades (' cluster_names{c} ' Cluster)'], ...
        'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'normal');

    legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], ...
           {'Young Pro', 'Young Anti', 'Old Pro', 'Old Anti'}, 'FontSize', 12, 'FontName', 'Helvetica');

    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

    xlim([-100, 500]);

    set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
    hold off;

    saveas(gcf, fullfile([' STIMULUS LOCKED Grand_Average_ERP_pro_vs_anti_' cluster_names{c} '.png']));
end



%% Plot for condition impact only (Pro vs Anti), collapsing across age
for c = 1:3
    figure('Color', 'white'); 
    hold on;

    % Calculate grand averages across age for Pro and Anti conditions
    pro_avg = (young_pro_avg{c} + old_pro_avg{c}) / 2;
    pro_sem = sqrt((young_pro_sem{c}.^2 + old_pro_sem{c}.^2) / 2);

    anti_avg = (young_anti_avg{c} + old_anti_avg{c}) / 2;
    anti_sem = sqrt((young_anti_sem{c}.^2 + old_anti_sem{c}.^2) / 2);

    % Plot Pro and Anti conditions with error bars
    h1 = shadedErrorBar(time, mean(pro_avg, 1), mean(pro_sem, 1), {'Color', [0.2, 0.6, 0.8], 'LineStyle', '-', 'LineWidth', 2});
    h2 = shadedErrorBar(time, mean(anti_avg, 1), mean(anti_sem, 1), {'Color', [0.8, 0.2, 0.2], 'LineStyle', '--', 'LineWidth', 2});

    xlabel('Time (ms)', 'FontSize', 14, 'FontName', 'Helvetica');
    ylabel('Amplitude (µV)', 'FontSize', 14, 'FontName', 'Helvetica');
    title(['Condition Impact: Pro vs Anti Saccades (' cluster_names{c} ' Cluster)'], ...
        'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'normal');

    legend([h1.mainLine, h2.mainLine], {'Pro', 'Anti'}, 'FontSize', 12, 'FontName', 'Helvetica');

    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

    xlim([-100, 500]);

    set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
    hold off;

    saveas(gcf, fullfile(['Condition_Impact_Pro_vs_Anti_' cluster_names{c} '.png']));
end

%% Plot for age impact only (Young vs Old), collapsing across conditions
for c = 1:3
    figure('Color', 'white'); 
    hold on;

    % Calculate grand averages across conditions for Young and Old groups
    young_avg = (young_pro_avg{c} + young_anti_avg{c}) / 2;
    young_sem = sqrt((young_pro_sem{c}.^2 + young_anti_sem{c}.^2) / 2);

    old_avg = (old_pro_avg{c} + old_anti_avg{c}) / 2;
    old_sem = sqrt((old_pro_sem{c}.^2 + old_anti_sem{c}.^2) / 2);

    % Plot Young and Old groups with error bars
    h1 = shadedErrorBar(time, mean(young_avg, 1), mean(young_sem, 1), {'Color', [0.3, 0.3, 0.9], 'LineStyle', '-', 'LineWidth', 2});
    h2 = shadedErrorBar(time, mean(old_avg, 1), mean(old_sem, 1), {'Color', [0.9, 0.3, 0.3], 'LineStyle', '--', 'LineWidth', 2});

    xlabel('Time (ms)', 'FontSize', 14, 'FontName', 'Helvetica');
    ylabel('Amplitude (µV)', 'FontSize', 14, 'FontName', 'Helvetica');
    title(['Age Impact: Young vs Old (' cluster_names{c} ' Cluster)'], ...
        'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'normal');

    legend([h1.mainLine, h2.mainLine], {'Young', 'Old'}, 'FontSize', 12, 'FontName', 'Helvetica');

    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

    xlim([-100, 500]);

    set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
    hold off;

    saveas(gcf, fullfile(['Age_Impact_Young_vs_Old_' cluster_names{c} '.png']));
end
