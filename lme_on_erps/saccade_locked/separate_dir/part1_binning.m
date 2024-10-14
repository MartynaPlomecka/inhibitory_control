% Single-Trial ERP Analysis Script with RT Binning and Averaging
eeglabPath = '\\psyger-stor02.d.uzh.ch\methlab\Students\Marius\toolboxes\eeglab2021.1';
HBNpath = '\\psyger-stor02.d.uzh.ch\methlab\HBN\EEG\Alpha analysis\FOOOF_HBN_final';
dataPath = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\dataT1_bas\test_baselined_clean_longsegmented_data_correcttrials_only_final_unfoldclean_opticat\' ;

if ~exist('pop_reref', 'file')
    addpath(eeglabPath);
    eeglab; close;
end

d = dir(dataPath);
d = d([d.isdir]);
d(startsWith({d.name}, '.')) = [];  

%%
for sj = 1:numel(d)
    try
        tic;
        disp(['Starting subject: ' num2str(sj)]);
        
        subjectPath = fullfile(d(sj).folder, d(sj).name);
        eegFile = dir(fullfile(subjectPath, '*sacclockedEEG.mat'));
        eegData = load(fullfile(eegFile.folder, eegFile.name));

        saccEEG = eegData.saccEEG;
        
        sacculus_events = saccEEG.event(2:2:end);
        age_group = sacculus_events(1).age; % 0 is young, 1 is old
        
        % Separate pro and anti trials, with left and right directions
        pro_left_indices = find(strcmp({sacculus_events.cond}, 'pro') & strcmp({sacculus_events.dir}, 'left'));
        pro_right_indices = find(strcmp({sacculus_events.cond}, 'pro') & strcmp({sacculus_events.dir}, 'right'));
        anti_left_indices = find(strcmp({sacculus_events.cond}, 'anti') & strcmp({sacculus_events.dir}, 'left'));
        anti_right_indices = find(strcmp({sacculus_events.cond}, 'anti') & strcmp({sacculus_events.dir}, 'right'));
        
        if isempty(pro_left_indices) || isempty(pro_right_indices) || isempty(anti_left_indices) || isempty(anti_right_indices)
            error('Not all pro/anti left/right trials found for subject %d', sj);
        end
        
        %  RT information
        pro_left_rt = [sacculus_events(pro_left_indices).rt];
        pro_right_rt = [sacculus_events(pro_right_indices).rt];
        anti_left_rt = [sacculus_events(anti_left_indices).rt];
        anti_right_rt = [sacculus_events(anti_right_indices).rt];
        
        % Bin trials and calculate mean RT for each bin
        [pro_left_binned, pro_left_bin_order, sorted_pro_left_rt, pro_left_bin_mean_rt] = bin_trials(pro_left_indices, pro_left_rt);
        [pro_right_binned, pro_right_bin_order, sorted_pro_right_rt, pro_right_bin_mean_rt] = bin_trials(pro_right_indices, pro_right_rt);
        [anti_left_binned, anti_left_bin_order, sorted_anti_left_rt, anti_left_bin_mean_rt] = bin_trials(anti_left_indices, anti_left_rt);
        [anti_right_binned, anti_right_bin_order, sorted_anti_right_rt, anti_right_bin_mean_rt] = bin_trials(anti_right_indices, anti_right_rt);
        
        disp('Binning for Pro Left trials:');
        log_binning_details(pro_left_binned, sorted_pro_left_rt, pro_left_bin_order);
        disp('Binning for Pro Right trials:');
        log_binning_details(pro_right_binned, sorted_pro_right_rt, pro_right_bin_order);
        disp('Binning for Anti Left trials:');
        log_binning_details(anti_left_binned, sorted_anti_left_rt, anti_left_bin_order);
        disp('Binning for Anti Right trials:');
        log_binning_details(anti_right_binned, sorted_anti_right_rt, anti_right_bin_order);
        
        % Average trials within bins
        pro_left_avg = cellfun(@(x) average_trials(saccEEG, x), pro_left_binned, 'UniformOutput', false);
        pro_right_avg = cellfun(@(x) average_trials(saccEEG, x), pro_right_binned, 'UniformOutput', false);
        anti_left_avg = cellfun(@(x) average_trials(saccEEG, x), anti_left_binned, 'UniformOutput', false);
        anti_right_avg = cellfun(@(x) average_trials(saccEEG, x), anti_right_binned, 'UniformOutput', false);
        
        % Stack averages into 3D arrays (channels x time points x bins)
        pro_left_erp = cat(3, pro_left_avg{:});
        pro_right_erp = cat(3, pro_right_avg{:});
        anti_left_erp = cat(3, anti_left_avg{:});
        anti_right_erp = cat(3, anti_right_avg{:});
        
        % Get time vector and channel locations
        time = saccEEG.times;  % in ms
        chanlocs = saccEEG.chanlocs;
        
        results_folder = fullfile(subjectPath, 'dir_ext_info_ERP_sacclock_binned_results');
        if ~exist(results_folder, 'dir')
            mkdir(results_folder);
        end
        save(fullfile(results_folder, 'ext_info_ERP_sacclock_binned_results.mat'), ...
            'pro_left_erp', 'pro_right_erp', 'anti_left_erp', 'anti_right_erp', 'time', 'chanlocs', 'age_group', ...
            'pro_left_rt', 'pro_right_rt', 'anti_left_rt', 'anti_right_rt', 'pro_left_binned', 'pro_right_binned', ...
            'anti_left_binned', 'anti_right_binned', 'pro_left_bin_mean_rt', 'pro_right_bin_mean_rt', ...
            'anti_left_bin_mean_rt', 'anti_right_bin_mean_rt', '-v7.3');
        
        disp(['Done subject ' num2str(sj) ' in ' num2str(toc) ' seconds']);

    catch ME
        warning(['Error with subject ' num2str(sj) ': ' ME.message']);
        disp(['Error details: ' getReport(ME)]);
    end
end

%% Functions

function [binned_indices, bin_order, sorted_rt, mean_rt] = bin_trials(indices, rt)

    % Sort RTs in ascending order and return corresponding indices
    [sorted_rt, sort_idx] = sort(rt);  
    sorted_indices = indices(sort_idx);  % Reorder the indices according to the sorted RTs
    n_trials = length(sorted_indices);
    bin_size = floor(n_trials / 3);  % 3 bins

    % Create bins with sorted indices and corresponding RTs
    binned_indices = {sorted_indices(1:bin_size), ...
                      sorted_indices(bin_size+1:2*bin_size), ...
                      sorted_indices(2*bin_size+1:end)};

    % Calculate mean RT for each bin
    mean_rt = [mean(sorted_rt(1:bin_size)), ...
               mean(sorted_rt(bin_size+1:2*bin_size)), ...
               mean(sorted_rt(2*bin_size+1:end))];

    % Bin order
    bin_order = {'fast', 'medium', 'slow'};
end

function log_binning_details(binned_indices, sorted_rt, bin_order)
    start_index = 1;
    for i = 1:length(binned_indices)
        disp(['Bin ' num2str(i) ' (' bin_order{i} '):']);
        disp(['  Number of trials: ' num2str(length(binned_indices{i}))]);
        disp('  Trials and corresponding RTs:');
        for j = 1:length(binned_indices{i})
            rt_value = sorted_rt(start_index);  
            trial_idx = binned_indices{i}(j);
            disp(['    Trial ' num2str(trial_idx) ': RT = ' num2str(rt_value)]);
            start_index = start_index + 1;
        end
    end
end

function avg_data = average_trials(saccEEG, trial_indices)
    % Average the data across trials specified by trial_indices
    avg_data = mean(saccEEG.data(:,:,trial_indices), 3);
end

