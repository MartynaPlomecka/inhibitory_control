%% RUN AFTER LME (PART 3)

coeffNames = lme_model.CoefficientNames;
disp(coeffNames);

cond_coeff = 'Condition_pro';
Mask = Results.(cond_coeff).Mask;

cond_idx = find(strcmp(lme_model.CoefficientNames, 'Condition_pro'));

% betas is a 105x301x6 matrix (channels x time x coefficients)
betas_to_use = betas(:, :, cond_idx);

time_windows = [-200, -100; -100 , -0];

time_indices = cell(size(time_windows, 1), 1);
for i = 1:size(time_windows, 1)
    time_indices{i} = find(time_points >= time_windows(i, 1) & time_points <= time_windows(i, 2));
end

% Averaging the betas within each time window
beta_averages = cell(size(time_windows, 1), 1);
for i = 1:length(time_indices)
    beta_averages{i} = mean(betas_to_use(:, time_indices{i}), 2);  % Avg across the selected time points
end

% Find significant electrodes in each time window (from Mask in Results)
significant_electrodes = cell(size(time_windows, 1), 1);
for i = 1:length(time_indices)
    significant_electrodes{i} = any(Mask(:, time_indices{i}), 2);
end


%% viz
figure('Color', 'w', 'InvertHardcopy', 'off', 'Position', [100, 100, 1200, 400]);  % Set a rectangular figure size
for i = 1:length(time_windows)
    subplot(1, 2, i);

    % Topoplot for the averaged betas in the current time window
    topoplot(beta_averages{i}, chanlocs, 'electrodes', 'on', 'style', 'both', 'conv', 'on', 'emarker2', {find(significant_electrodes{i}), 'o', 'w', 3, 0.7});

    colormap(mycolormap);

    c = colorbar;
    c.Position = [c.Position(1) + 0.05, c.Position(2), c.Position(3) * 0.8, c.Position(4)];
    clim([-0.1 0.1]);

    title(['Time: ' num2str(time_windows(i, 1)) '-' num2str(time_windows(i, 2)) ' ms'], 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Helvetica');
end

sgtitle('Condition Coefficient', 'FontSize', 18, 'FontWeight', 'normal', 'FontName', 'Helvetica');
