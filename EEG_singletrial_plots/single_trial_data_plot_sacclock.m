cd \\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_new_baselineremov_unfoldclean\AB2
load('AB2_sacclockedEEG.mat')
load('AB2_stimlockedEEG.mat')

%electrodes to average
electrodes_to_average = [58, 62, 69];

% Choose 
conditionType = 'PRO'; % Or 'ANTI'

[Cond, Cond_rt] = selectCondition(conditionType, saccEEG);

%mean EEG data for selected electrodes and trials
dat_Cond = squeeze(mean(saccEEG.data(electrodes_to_average, :, Cond), 1));
%sorting data based on RTs
[sorted, sortInd] = sort(Cond_rt);
dat_Cond_sorted = dat_Cond(:, sortInd);
% index for time zero
time_zero_index = find(saccEEG.times == 0);

%% Plotting
figure('Color', [1 1 1]); 
colormap(customcolormap_preset('red-white-blue')); 

% Subplot for unsorted data
subplot(121)
imagesc(dat_Cond')
clim([-40 40]); 
colorbar; 
xticks(1:100:length(dat_Cond'))
xticklabels(saccEEG.times(1:100:end))
title([conditionType ' condition: unsorted'], 'FontSize', 14, 'FontWeight', 'normal')

% Subplot for sorted data with StimOnset and RTs
subplot(122)
imagesc(dat_Cond_sorted')
clim([-40 40]); 
colorbar; 
xticks(1:100:length(dat_Cond_sorted'))
xticklabels(saccEEG.times(1:100:end))
title([conditionType ' condition: sorted with StimOnset'], 'FontSize', 14, 'FontWeight', 'normal')
hold on;

%  vertical line at the saccade onset
line([time_zero_index, time_zero_index], ylim, 'Color', 'k', 'LineWidth', 2);

% Add reaction times 
dotSize = 5; 
for i = 1:length(Cond)
    x_coord = time_zero_index - (sorted(i) / 2); % Adjusting based on sampling rate (500 Hz)
    scatter(x_coord, i, dotSize, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
end

hold off;

function [Cond, Cond_rt] = selectCondition(conditionType, saccEEG)
    % Select condition (Pro or Anti) and return corresponding epochs and RTs
    switch conditionType
        case 'PRO'
            Cond = [saccEEG.event(find(strcmp({saccEEG.event.type}, '22') | strcmp({saccEEG.event.type}, '21'))).epoch];
            Cond_rt = [saccEEG.event(find(strcmp({saccEEG.event.type}, '22') | strcmp({saccEEG.event.type}, '21'))).rt];
        case 'ANTI'
            Cond = [saccEEG.event(find(strcmp({saccEEG.event.type}, '24') | strcmp({saccEEG.event.type}, '23'))).epoch];
            Cond_rt = [saccEEG.event(find(strcmp({saccEEG.event.type}, '24') | strcmp({saccEEG.event.type}, '23'))).rt];
        otherwise
            error('Invalid condition type. Choose either ''PRO'' or ''ANTI''.');
    end
end
