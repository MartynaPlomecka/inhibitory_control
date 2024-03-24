mainDir = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_new_baselineremov_unfoldclean';

folderInfo = dir(mainDir);
folderNames = {folderInfo([folderInfo.isdir]).name};
subjectDirs = folderNames(~ismember(folderNames, {'.', '..'})); % Exclude '.' and '..'

combined_dat_Cond_young = [];
combined_Cond_rt_young = [];
combined_dat_Cond_old = [];
combined_Cond_rt_old = [];

conditionType = 'ANTI'; %  'PRO' or 'ANTI'

for i = 1:10 %length(subjectDirs)
    subjectDir = fullfile(mainDir, subjectDirs{i});
    eegFile = fullfile(subjectDir, [subjectDirs{i}, '_sacclockedEEG.mat']);
    
    if exist(eegFile, 'file')
        load(eegFile); 
        electrodes_to_average = [8, 13, 14];

        isOld = saccEEG.event(1).age; % 0 for young, 1 for old
        [Cond, Cond_rt] = selectCondition(conditionType, saccEEG);
        dat_Cond = squeeze(mean(saccEEG.data(electrodes_to_average, :, Cond), 1));

        if isOld
            combined_dat_Cond_old = cat(2, combined_dat_Cond_old, dat_Cond);
            combined_Cond_rt_old = [combined_Cond_rt_old; Cond_rt'];
        else
            combined_dat_Cond_young = cat(2, combined_dat_Cond_young, dat_Cond);
            combined_Cond_rt_young = [combined_Cond_rt_young; Cond_rt']; 
        end
    else
        warning(['File does not exist: ', eegFile]);
    end
end

% Plot for young subjects
plotData(combined_dat_Cond_young, combined_Cond_rt_young, saccEEG.times, 'Young', 'red-white-blue', conditionType);

% Plot for old subjects
plotData(combined_dat_Cond_old, combined_Cond_rt_old, saccEEG.times, 'Old', 'red-white-blue', conditionType);

% Function definitions
function plotData(data, rt, times, ageGroup, colormapName, conditionType)
    [rt_sorted, sortIndex] = sort(rt);
    data_sorted = data(:, sortIndex);

    figure('Color', [1 1 1]); 
    colormap(customcolormap_preset(colormapName)); 

    % Subplot for unsorted data
    subplot(121)
    imagesc(data')
    clim([-30 30]); 
    colorbar; 
    xticks(1:100:length(data))
    xticklabels(times(1:100:end))
    title([ageGroup ' ' conditionType ' condition: unsorted'], 'FontSize', 14, 'FontWeight', 'normal')

    % Subplot for sorted data
    subplot(122)
    imagesc(data_sorted')
    clim([-30 30]); 
    colorbar; 
    xticks(1:100:length(data_sorted))
    xticklabels(times(1:100:end))
    title([ageGroup ' ' conditionType ' condition: sorted'], 'FontSize', 14, 'FontWeight', 'normal')
    hold on;

    % Vertical line at the saccade onset
    time_zero_index = find(times == 0);
    line([time_zero_index, time_zero_index], ylim, 'Color', 'k', 'LineWidth', 2);

    % Add reaction times
    dotSize = 5; 
    for i = 1:length(rt_sorted)
        x_coord = time_zero_index - (rt_sorted(i) / 2);
        scatter(x_coord, i, dotSize, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    end

    hold off;
end

function [Cond, Cond_rt] = selectCondition(conditionType, saccEEG)
    switch conditionType
        case 'PRO'
            Cond = [saccEEG.event(strcmp({saccEEG.event.type}, '22') | strcmp({saccEEG.event.type}, '21')).epoch];
            Cond_rt = [saccEEG.event(strcmp({saccEEG.event.type}, '22') | strcmp({saccEEG.event.type}, '21')).rt];
        case 'ANTI'
            Cond = [saccEEG.event(strcmp({saccEEG.event.type}, '24') | strcmp({saccEEG.event.type}, '23')).epoch];
            Cond_rt = [saccEEG.event(strcmp({saccEEG.event.type}, '24') | strcmp({saccEEG.event.type}, '23')).rt];
        otherwise
            error('Invalid condition type. Choose either ''PRO'' or ''ANTI''.');
    end
end
