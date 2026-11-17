%% Mean Spectrum + MTRasym for 5 Animals (each averaged from 3 repeats)
% Data: TestData_Long_Ctrl_HPD_CESTBase.csv
% Structure: subplot(5,2) - 5 rows (animals), col1: Mean+SD spectrum, col2: MTRasym mean

% Load data
data = csvread('TestData_Long_Ctrl_HFD_CESTBase.csv', 1, 0);  % 33 rows x 15 cols
freq_original = linspace(-5, 5, 33);
xaxis = linspace(-5, 5, 4000);

% Define animals and their repeat columns (1-based indices)
% animals = {'802', '803', '804', '805', '806'};
animals = {'816', '818','819'};

col_groups = {[1,2,3], [4,5,6],[7,8,9]};  % Column indices for each animal's 3 repeats
% col_groups = {[1,2,3], [4,5,6], [7,8,9], [10,11,12], [13,14,15]};  % Column indices for each animal's 3 repeats

n_animals = 3;

% Create figure with 5 rows x 2 columns subplots
figure('Position', [50 50 1400 2000]);

for a = 1:n_animals
    cols = col_groups{a};
    n_repeats = 3;
    
    % Storage for B0-corrected spectra for this animal's repeats
    z_corrected_repeats = zeros(n_repeats, length(xaxis));
    
    % B0 correct each repeat
    for r = 1:n_repeats
        z_spectrum = data(:, cols(r));
        
        % Spline to 4000 points
        freq_spline = linspace(freq_original(1), freq_original(end), 4000);
        z_spline = spline(freq_original, z_spectrum, freq_spline);
        
        % Find B0 minimum
        [~, min_idx] = min(z_spline);
        b0_shift = freq_spline(min_idx);
        
        % Shift frequency axis
        freq_shifted = freq_spline - b0_shift;
        
        % Evaluate on common grid
        pp1 = spline(freq_shifted, z_spline);
        z_corrected_repeats(r, :) = ppval(pp1, xaxis);
    end
    
    % Average across repeats for this animal
    mean_spectrum = mean(z_corrected_repeats, 1);
    std_spectrum = std(z_corrected_repeats, 0, 1);  % Population std across repeats
    
    % MTRasym spectrum (mean only, no shade)
    mid_idx = round(length(xaxis)/2);
    freq_positive = xaxis(mid_idx:end);
    MTRasym_spectrum = zeros(size(freq_positive));
    
    for j = 1:length(freq_positive)
        freq_pos = freq_positive(j);
        freq_neg = -freq_pos;
        idx_pos = find(xaxis >= freq_pos, 1);
        idx_neg = find(xaxis <= freq_neg, 1, 'last');
        if ~isempty(idx_pos) && ~isempty(idx_neg)
            MTRasym_spectrum(j) = mean_spectrum(idx_neg) - mean_spectrum(idx_pos);
        end
    end
    
    % MTRasym at 1.28 ppm
    idx_128 = find(freq_positive >= 1.28, 1);
    mtrasym_128 = MTRasym_spectrum(idx_128);
    fprintf('%s: MTRasym(1.28 ppm) = %.6f\n', animals{a}, mtrasym_128);
    
    %% Subplot 1: Mean spectrum + shaded SD (across 3 repeats)
    subplot(n_animals, 2, (a-1)*2 + 1);
    hold on;
    fill([xaxis, fliplr(xaxis)], [mean_spectrum + std_spectrum, fliplr(mean_spectrum - std_spectrum)], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xaxis, mean_spectrum, 'b-', 'LineWidth', 2.5);
    plot([1.28 1.28], [min(mean_spectrum) max(mean_spectrum)], 'r--', 'LineWidth', 1.5);
%     plot([-1.28 -1.28], [min(mean_spectrum) max(mean_spectrum)], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Frequency offset (ppm)');
    ylabel('Z/Z_0');
    title(sprintf('%s (n=3 repeats): Mean CEST Spectrum', animals{a}));
    grid on;
    xlim([-5 5]);
    ylim([0 1]);  % Adjust based on data
    set(gca, 'XDir', 'reverse');
    hold off;
    
    %% Subplot 2: MTRasym mean (NO shading)
    subplot(n_animals, 2, (a-1)*2 + 2);
    hold on;
    plot(freq_positive, MTRasym_spectrum, 'r-', 'LineWidth', 2.5);
    plot([1.28 1.28], [min(MTRasym_spectrum) max(MTRasym_spectrum)], 'k--', 'LineWidth', 1.5);
    plot([0 5], [0 0], 'k:', 'LineWidth', 1);
    xlabel('Frequency offset (ppm)');
    ylabel('MTRasym');
    title(sprintf('%s: Asymmetric Spectrum', animals{a}));
    grid on;
    xlim([0 5]);
    ylim([-0.08 0.03]);  % Adjust as needed
    set(gca, 'XDir', 'reverse');
    hold off;
end

