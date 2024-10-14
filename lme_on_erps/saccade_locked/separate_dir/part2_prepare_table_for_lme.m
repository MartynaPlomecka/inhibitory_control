clc;
clear;

dataPath = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_newest\THETAproject\dataT1_bas\test_baselined_clean_longsegmented_data_correcttrials_only_final_unfoldclean_opticat\' ;
speed_bins = {'fast', 'medium', 'slow'};

d = dir(dataPath);
d = d([d.isdir]);
d(startsWith({d.name}, '.')) = [];

all_subject_tables = table();  
for sj = 1:numel(d)
    try
        subjectPath = fullfile(d(sj).folder, d(sj).name);
        resultsFile = fullfile(subjectPath, 'dir_ext_info_ERP_sacclock_binned_results', 'ext_info_ERP_sacclock_binned_results.mat');
        
        load(resultsFile);
        subject_id = d(sj).name;
        conditions = {'pro', 'anti'};  
        directions = {'left', 'right'};  
        subject_table = table();  

        for cond_idx = 1:2
            cond_name = conditions{cond_idx};
            for dir_idx = 1:2
                dir_name = directions{dir_idx};
                cond_dir_erp = eval([cond_name '_' dir_name '_erp']);  
                cond_dir_mean_rt = eval([cond_name '_' dir_name '_bin_mean_rt']);  
                age_group = eval('age_group');  %  (0 = young, 1 = old)

                % Loop through bins (fast, medium, slow)
                for bin_idx = 1:3
                    speed_bin = speed_bins{bin_idx};  % speed bin (fast, medium, slow)
                    mean_rt = cond_dir_mean_rt(bin_idx);  %  mean RT for the bin
                    eeg_data = cond_dir_erp(:,:,bin_idx);  % ERP data for the current bin (channels x time points)

                    new_row = table({subject_id}, age_group, {cond_name}, {dir_name}, {speed_bin}, mean_rt, {eeg_data}, ...
                        'VariableNames', {'SubjectID', 'AgeGroup', 'Condition', 'Direction', 'SpeedBin', 'MeanRT', 'EEGData'});

                    subject_table = [subject_table; new_row];
                end
            end
        end
        
        all_subject_tables = [all_subject_tables; subject_table];
        
        disp(['Processed subject ' subject_id]);

    catch ME
        warning(['Error processing subject ' d(sj).name ': ' ME.message']);
    end
end

save('all_subjects_sacclocked_ERP_lme_table_separated_direction.mat', 'all_subject_tables');
