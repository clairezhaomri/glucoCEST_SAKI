%% Average 6 spectra and plot mean asymmetric spectrum
% Load data
data = csvread('TestData_CtrlAcuteSCCESTBase.csv', 1, 0);
animals = {'11044', '800', '801', '880', '945'};
n_animals = 5;
freq_original = linspace(-5, 5, 33);

% Define common x-axis
xaxis = linspace(-5, 5, 4000);

% Storage for B0-corrected spectra
z_corrected_all = zeros(n_animals, length(xaxis));

% B0 correct each animal and store
for i = 1:n_animals
    z_spectrum = data(:, i);
    
    % Spline to 4000 points
    freq_spline = linspace(freq_original(1), freq_original(end), 4000);
    z_spline = spline(freq_original, z_spectrum, freq_spline);
    
    % Find B0 minimum
    [~, min_idx] = min(z_spline);
    b0_shift = freq_spline(min_idx);
    
    % Shift frequency axis
    freq_shifted = freq_spline - b0_shift;
    
    % Create spline and evaluate on common grid
    pp1 = spline(freq_shifted, z_spline);
    z_corrected_all(i, :) = ppval(pp1, xaxis);
end

% Calculate mean spectrum
mean_spectrum = mean(z_corrected_all, 1);
std_spectrum = std(z_corrected_all, 0, 1);

% Calculate MTRasym spectrum
% MTRasym(ω) = Z(-ω) - Z(+ω)
% Need to find corresponding negative frequencies for each positive frequency
n_points = length(xaxis);
mid_idx = round(n_points/2);

% Create symmetric frequency axis (only positive side)
freq_positive = xaxis(mid_idx:end);
MTRasym_spectrum = zeros(size(freq_positive));

for j = 1:length(freq_positive)
    freq_pos = freq_positive(j);
    freq_neg = -freq_pos;
    
    % Find indices
    idx_pos = find(xaxis >= freq_pos, 1);
    idx_neg = find(xaxis <= freq_neg, 1, 'last');
    
    if ~isempty(idx_pos) && ~isempty(idx_neg)
        MTRasym_spectrum(j) = mean_spectrum(idx_neg) - mean_spectrum(idx_pos);
    end
end

%% Plot 1: Mean spectrum with ±1 std
figure('Position', [100 100 1200 400]);

subplot(1,2,1);
hold on;
% Plot mean with shaded std
fill([xaxis, fliplr(xaxis)], ...
     [mean_spectrum + std_spectrum, fliplr(mean_spectrum - std_spectrum)], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(xaxis, mean_spectrum, 'b-', 'LineWidth', 2.5);
plot([1.28 1.28], [0.55 0.75], 'r--', 'LineWidth', 1.5);
plot([-1.28 -1.28], [0.55 0.75], 'r--', 'LineWidth', 1.5);
xlabel('Frequency offset (ppm)');
ylabel('Normalized Signal (Z/Z_0)');
title('Mean CEST Spectrum (n=6)');
grid on;
xlim([-5 5]);
ylim([0 1]);
hold off;

%% Plot 2: MTRasym spectrum
subplot(1,2,2);
hold on;
plot(freq_positive, MTRasym_spectrum, 'r-', 'LineWidth', 2.5);
plot([1.28 1.28], [min(MTRasym_spectrum) max(MTRasym_spectrum)], ...
     'k--', 'LineWidth', 1.5);
plot([0 5], [0 0], 'k:', 'LineWidth', 1);
xlabel('Frequency offset (ppm)');
ylabel('MTRasym');
title('Asymmetric Spectrum');
grid on;
xlim([0 5]);
ylim([-0.08 0.03]);
hold off;

% Print MTRasym value at 1.28 ppm
idx_128 = find(freq_positive >= 1.28, 1);
fprintf('\nMTRasym at 1.28 ppm = %.6f\n', MTRasym_spectrum(idx_128));

%% Plot 3: All individual spectra + mean
figure('Position', [100 100 1000 600]);
colors = lines(n_animals);
hold on;
for i = 1:n_animals
    plot(xaxis, z_corrected_all(i,:), 'Color', [colors(i,:) 0.3], ...
         'LineWidth', 1, 'DisplayName', animals{i});
end
plot(xaxis, mean_spectrum, 'k-', 'LineWidth', 3, 'DisplayName', 'Mean');
plot([1.28 1.28], [0.55 0.75], 'r--', 'LineWidth', 3, 'DisplayName', '±1.28 ppm');
plot([-1.28 -1.28], [0.55 0.75], 'r--', 'LineWidth',3, 'HandleVisibility', 'off');
xlabel('Frequency offset (ppm)');
ylabel('Normalized Signal (Z/Z_0)');
title('All Individual Spectra + Mean');
grid on;
xlim([-5 5]);
ylim([0 1]);
hold off;
