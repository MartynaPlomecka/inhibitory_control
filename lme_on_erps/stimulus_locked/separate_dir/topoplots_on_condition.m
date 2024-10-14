coeff = 'Condition_pro';

% Get Obs (== betas?) and the significance mask
Obs = Results.(coeff).Obs; 
Mask = Results.(coeff).Mask;

% Define time windows (in ms)
time_windows = [200, 250; 300, 350; 450, 500];

time_indices = cell(size(time_windows, 1), 1);
for i = 1:size(time_windows, 1)
    time_indices{i} = find(time_points >= time_windows(i, 1) & time_points <= time_windows(i, 2));
end

% Averaging the betas within each time window
beta_averages = cell(size(time_windows, 1), 1);
for i = 1:length(time_indices)
    beta_averages{i} = mean(Obs(:, time_indices{i}), 2);  % Average across the selected time points
end

% Find significant electrodes in each time window
significant_electrodes = cell(size(time_windows, 1), 1);
for i = 1:length(time_indices)
    significant_electrodes{i} = any(Mask(:, time_indices{i}), 2);  % Check if any time point in the window is significant
end

% Create the topoplots for each time window
figure('Color', 'w', 'InvertHardcopy', 'off', 'Position', [100, 100, 1200, 400]);  % Set a rectangular figure size
for i = 1:length(time_windows)
    subplot(1, 3, i);
    
    % Topoplot for the averaged betas in the current time window
    topoplot(beta_averages{i}, chanlocs, 'electrodes', 'on', 'style', 'both', 'conv', 'on', 'emarker2', {find(significant_electrodes{i}), 'o', 'w', 3, 0.7});  
    
    colormap(mycolormap);
    
    c = colorbar;
    c.Position = [c.Position(1) + 0.05, c.Position(2), c.Position(3) * 0.8, c.Position(4)];  
    clim([-4, 4]);  
    
    title(['Time: ' num2str(time_windows(i, 1)) '-' num2str(time_windows(i, 2)) ' ms'], 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Helvetica');
end

sgtitle('Condition coefficient', 'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica');  