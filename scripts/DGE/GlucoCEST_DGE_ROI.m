function [AnimalID, DeltaS_Result] = GlucoCEST_DGE_ROI(pathBase, pathInfu, SelectPpm, ROI_name)
% DGE-style glucoCEST with PIXELWISE B0 correction + FWHM spillover correction
% Modified Feb 2026 - Complete version with MTRasym analysis
% 
% Inputs:
%   pathBase: path to baseline CEST folder
%   pathInfu: path to infusion CEST folder
%   SelectPpm: single frequency (e.g., 1.28) or range [min, max]
%   ROI_name: 'ROI_ctx' or 'ROI_medulla'

%% Load baseline data
CESTbase = MRIDataLoader.LoadFrom2DSeq([pathBase filesep 'pdata' filesep '1' filesep '2dseq']);
AnimalID = CESTbase.header.patientID;

if isfile([pathBase filesep 'CESTCorrectData.mat'])
    fprintf('  Loading motion-corrected base data from CESTCorrectData.mat\n');
    temp = load([pathBase filesep 'CESTCorrectData.mat'], 'CESTCorrectData');
    
    if isa(temp.CESTCorrectData, 'MRIData')
        corrected_data = temp.CESTCorrectData.vol;
    else
        corrected_data = temp.CESTCorrectData;
    end
    
    if isequal(size(CESTbase.vol), size(corrected_data))
        CESTbase.vol = corrected_data;
        fprintf('    Successfully replaced with motion-corrected data\n');
    else
        warning('Dimension mismatch! Using raw data');
    end
end

CESTbase.Permute('echo','read','line');

if strcmp(CESTbase.headerHistory{1}.method.M0Option,'M0_Last')
    S0_base = squeeze(CESTbase.vol(CESTbase.Size('echo'),:,:));
    CESTbase.vol(CESTbase.Size('echo'),:,:) = [];
elseif strcmp(CESTbase.headerHistory{1}.method.M0Option,'M0_First')
    S0_base = squeeze(CESTbase.vol(1,:,:));
    CESTbase.vol(1,:,:) = [];
end

%% Load infusion data
CESTinfu = MRIDataLoader.LoadFrom2DSeq([pathInfu filesep 'pdata' filesep '1' filesep '2dseq']);

if isfile([pathInfu filesep 'CESTCorrectData.mat'])
    fprintf('  Loading motion-corrected infu data from CESTCorrectData.mat\n');
    temp = load([pathInfu filesep 'CESTCorrectData.mat'], 'CESTCorrectData');
    
    if isa(temp.CESTCorrectData, 'MRIData')
        corrected_data = temp.CESTCorrectData.vol;
    else
        corrected_data = temp.CESTCorrectData;
    end
    
    if isequal(size(CESTinfu.vol), size(corrected_data))
        CESTinfu.vol = corrected_data;
        fprintf('    Successfully replaced with motion-corrected data\n');
    else
        warning('Dimension mismatch! Using raw data');
    end
end

CESTinfu.Permute('echo','read','line');

if strcmp(CESTinfu.headerHistory{1}.method.M0Option,'M0_Last')
    S0_infu = squeeze(CESTinfu.vol(CESTinfu.Size('echo'),:,:));
    CESTinfu.vol(CESTinfu.Size('echo'),:,:) = [];
elseif strcmp(CESTinfu.headerHistory{1}.method.M0Option,'M0_First')
    S0_infu = squeeze(CESTinfu.vol(1,:,:));
    CESTinfu.vol(1,:,:) = [];
end

%% Load ROI mask
if isfile([pathBase filesep ROI_name '.mat'])
    roi_data = load([pathBase filesep ROI_name '.mat']);
    fprintf('  Loaded ROI mask from base folder: %s.mat\n', ROI_name);
elseif isfile([pathInfu filesep ROI_name '.mat'])
    roi_data = load([pathInfu filesep ROI_name '.mat']);
    fprintf('  Loaded ROI mask from infu folder: %s.mat\n', ROI_name);
else
    error('ROI mask %s.mat not found!', ROI_name);
end

field_names = fieldnames(roi_data);
if ismember('ROI', field_names)
    ROI_mask = roi_data.ROI;
elseif length(field_names) == 1
    ROI_mask = roi_data.(field_names{1});
else
    error('Cannot determine ROI mask variable');
end

if ~islogical(ROI_mask)
    ROI_mask = logical(ROI_mask);
end

fprintf('  ROI contains %d voxels\n', sum(ROI_mask(:)));

%% Get offsets and validate
MT_Offsets = CESTbase.headerHistory{1}.method.MT_Offsets_NoM0(:).';
nOffsets = length(MT_Offsets);

if isempty(nOffsets) || nOffsets < 5
    error('Invalid number of offsets: %d', nOffsets);
end

if size(CESTbase.vol, 1) ~= nOffsets || size(CESTinfu.vol, 1) ~= nOffsets
    error('Data size mismatch!');
end

fprintf('  Found %d frequency offsets\n', nOffsets);

%% PIXELWISE B0 CORRECTION (like your method)
fprintf('  Performing PIXELWISE B0 correction for baseline...\n');

[roi_x, roi_y] = find(ROI_mask);
nROIpixels = length(roi_x);
fprintf('    Processing %d pixels in ROI...\n', nROIpixels);

CESTbase_B0corrected = zeros(size(CESTbase.vol));

for p = 1:nROIpixels
    x = roi_x(p);
    y = roi_y(p);
    
    pixel_spectrum = squeeze(CESTbase.vol(:, x, y)) ./ S0_base(x, y);
    
    if any(~isfinite(pixel_spectrum))
        continue;
    end
    
    [~, min_idx] = min(pixel_spectrum);
    w0_pixel = MT_Offsets(min_idx);
    offsets_corrected = MT_Offsets - w0_pixel;
    
    [offsets_sorted, sort_idx] = sort(offsets_corrected);
    spectrum_sorted = pixel_spectrum(sort_idx);
    
    pp_pixel = spline(offsets_sorted, spectrum_sorted);
    CESTbase_B0corrected(:, x, y) = ppval(pp_pixel, MT_Offsets);
end

fprintf('  Performing PIXELWISE B0 correction for infusion...\n');
CESTinfu_B0corrected = zeros(size(CESTinfu.vol));

for p = 1:nROIpixels
    x = roi_x(p);
    y = roi_y(p);
    
    pixel_spectrum = squeeze(CESTinfu.vol(:, x, y)) ./ S0_infu(x, y);
    
    if any(~isfinite(pixel_spectrum))
        continue;
    end
    
    [~, min_idx] = min(pixel_spectrum);
    w0_pixel = MT_Offsets(min_idx);
    offsets_corrected = MT_Offsets - w0_pixel;
    
    [offsets_sorted, sort_idx] = sort(offsets_corrected);
    spectrum_sorted = pixel_spectrum(sort_idx);
    
    pp_pixel = spline(offsets_sorted, spectrum_sorted);
    CESTinfu_B0corrected(:, x, y) = ppval(pp_pixel, MT_Offsets);
end

%% Average B0-corrected pixels within ROI
fprintf('  Averaging B0-corrected spectra within ROI...\n');

Z_base = zeros(1, nOffsets);
Z_infu = zeros(1, nOffsets);

for i = 1:nOffsets
    slice_base = squeeze(CESTbase_B0corrected(i, :, :));
    slice_infu = squeeze(CESTinfu_B0corrected(i, :, :));
    
    Z_base(i) = mean(slice_base(ROI_mask), 'omitnan');
    Z_infu(i) = mean(slice_infu(ROI_mask), 'omitnan');
end

%% Fine-tune B0 on averaged spectrum
[~, xIndex] = sort(MT_Offsets);
CEST_raw_base = [MT_Offsets(xIndex); Z_base(xIndex)];

[~, min_idx] = min(CEST_raw_base(2,:));
w0_base_residual = CEST_raw_base(1, min_idx);
CEST_corrected_base = [CEST_raw_base(1,:) - w0_base_residual; CEST_raw_base(2,:)];

x_spline = linspace(CEST_corrected_base(1,1), CEST_corrected_base(1,end), 4000);
y_spline_base = spline(CEST_corrected_base(1,:), CEST_corrected_base(2,:), x_spline);

[~, min_idx_spline] = min(y_spline_base);
w0_spline_base = x_spline(min_idx_spline);
x_final_base = x_spline - w0_spline_base;
pp_base = spline(x_final_base, y_spline_base);

fprintf('  Base residual B0 shift: %.3f ppm\n', w0_spline_base/400);

CEST_raw_infu = [MT_Offsets(xIndex); Z_infu(xIndex)];
[~, min_idx] = min(CEST_raw_infu(2,:));
w0_infu_residual = CEST_raw_infu(1, min_idx);
CEST_corrected_infu = [CEST_raw_infu(1,:) - w0_infu_residual; CEST_raw_infu(2,:)];

y_spline_infu = spline(CEST_corrected_infu(1,:), CEST_corrected_infu(2,:), x_spline);
[~, min_idx_spline] = min(y_spline_infu);
w0_spline_infu = x_spline(min_idx_spline);
x_final_infu = x_spline - w0_spline_infu;
pp_infu = spline(x_final_infu, y_spline_infu);

fprintf('  Infu residual B0 shift: %.3f ppm\n\n', w0_spline_infu/400);

%% FWHM correction with fallback
fprintf('  Calculating FWHM for spillover correction...\n');

[FWHM_base, min_base, success_base] = calculate_water_FWHM(pp_base, 400);
[FWHM_infu, min_infu, success_infu] = calculate_water_FWHM(pp_infu, 400);

fprintf('  Water peak: base_min=%.4f, infu_min=%.4f\n', min_base, min_infu);

if success_base && success_infu
    fprintf('  Water FWHM: base=%.3f ppm, infu=%.3f ppm\n', FWHM_base, FWHM_infu);
    spillover_ratio = FWHM_base / FWHM_infu;
    fprintf('  Spillover correction ratio: %.3f\n', spillover_ratio);
    
    x_common_Hz = linspace(-2000, 2000, 4000);
    y_infu_raw = ppval(pp_infu, x_common_Hz);
    y_infu_corrected = (y_infu_raw - min_infu) * spillover_ratio + min_infu;
    pp_infu_corrected = spline(x_common_Hz, y_infu_corrected);
    correction_applied = true;
    
elseif ~success_base && success_infu
    fprintf('  WARNING: Base FWHM failed, using water depth ratio\n');
    fprintf('  Water FWHM: base=FAILED, infu=%.3f ppm\n', FWHM_infu);
    
    depth_ratio = (1 - min_base) / (1 - min_infu);
    spillover_ratio = depth_ratio;
    fprintf('  Using water depth ratio: %.3f\n', spillover_ratio);
    
    x_common_Hz = linspace(-2000, 2000, 4000);
    y_infu_raw = ppval(pp_infu, x_common_Hz);
    y_infu_corrected = (y_infu_raw - min_infu) * spillover_ratio + min_infu;
    pp_infu_corrected = spline(x_common_Hz, y_infu_corrected);
    correction_applied = true;
    FWHM_base = NaN;
    
else
    fprintf('  WARNING: FWHM correction failed, skipping\n');
    pp_infu_corrected = pp_infu;
    spillover_ratio = 1.0;
    correction_applied = false;
    if ~success_base; FWHM_base = NaN; end
    if ~success_infu; FWHM_infu = NaN; end
end

%% Calculate ΔS
if length(SelectPpm) == 1
    SelectFreq_Hz = SelectPpm * 400;
    S_base_at_freq = ppval(pp_base, SelectFreq_Hz);
    S_infu_at_freq = ppval(pp_infu_corrected, SelectFreq_Hz);
    
    DeltaS_Result.type = 'single_offset';
    DeltaS_Result.offset_ppm = SelectPpm;
    DeltaS_Result.S_base = S_base_at_freq;
    DeltaS_Result.S_infu_corrected = S_infu_at_freq;
    DeltaS_Result.S_infu_uncorrected = ppval(pp_infu, SelectFreq_Hz);
    DeltaS_Result.DeltaS_percent = (S_base_at_freq - S_infu_at_freq) * 100;
    
    fprintf('\n=== Results (%.2f ppm) ===\n', SelectPpm);
    fprintf('  S_base = %.4f\n', S_base_at_freq);
    fprintf('  S_infu (uncorrected) = %.4f\n', DeltaS_Result.S_infu_uncorrected);
    fprintf('  S_infu (corrected) = %.4f\n', S_infu_at_freq);
    fprintf('  ΔS (corrected) = %.2f%%\n', DeltaS_Result.DeltaS_percent);
    
else
    freq_min = SelectPpm(1) * 400;
    freq_max = SelectPpm(2) * 400;
    x_integration = linspace(freq_min, freq_max, 100);
    
    y_base_band = ppval(pp_base, x_integration);
    y_infu_band = ppval(pp_infu_corrected, x_integration);
    
    S_base_mean = trapz(x_integration, y_base_band) / (freq_max - freq_min);
    S_infu_mean = trapz(x_integration, y_infu_band) / (freq_max - freq_min);
    
    DeltaS_Result.type = 'band_integration';
    DeltaS_Result.ppm_range = SelectPpm;
    DeltaS_Result.S_base_mean = S_base_mean;
    DeltaS_Result.S_infu_mean = S_infu_mean;
    DeltaS_Result.DeltaS_percent = (S_base_mean - S_infu_mean) * 100;
    
    fprintf('\n=== Results (%.1f-%.1f ppm band) ===\n', SelectPpm(1), SelectPpm(2));
    fprintf('  ΔS (corrected) = %.2f%%\n', DeltaS_Result.DeltaS_percent);
end

%% Save metadata
DeltaS_Result.FWHM_base = FWHM_base;
DeltaS_Result.FWHM_infu = FWHM_infu;
DeltaS_Result.spillover_ratio = spillover_ratio;
DeltaS_Result.correction_applied = correction_applied;
DeltaS_Result.water_min_base = min_base;
DeltaS_Result.water_min_infu = min_infu;
DeltaS_Result.Z_base_raw = Z_base;
DeltaS_Result.Z_infu_raw = Z_infu;
DeltaS_Result.offsets_raw = MT_Offsets;
DeltaS_Result.pp_base = pp_base;
DeltaS_Result.pp_infu_corrected = pp_infu_corrected;
DeltaS_Result.ROI_name = ROI_name;
DeltaS_Result.ROI_mask = ROI_mask;
DeltaS_Result.AnimalID = AnimalID;

%% Visualization (4 subplots) - FINAL VERSION with fixed legends
figure('Position', [50 100 1600 350]);

%% Subplot 1: Raw data points + Spline fits
subplot(1,4,1);

% Convert Hz to ppm for raw data
[offsets_sorted_Hz, sort_idx] = sort(MT_Offsets);
offsets_sorted_ppm = offsets_sorted_Hz / 400;

% Plot raw measurement points
h1 = plot(offsets_sorted_ppm, Z_base(sort_idx), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
h2 = plot(offsets_sorted_ppm, Z_infu(sort_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% Plot spline fits
x_plot_ppm = linspace(-4, 4, 1000);
x_plot_Hz = x_plot_ppm * 400;
y_base_plot = ppval(pp_base, x_plot_Hz);
y_infu_plot = ppval(pp_infu, x_plot_Hz);

h3 = plot(x_plot_ppm, y_base_plot, 'b-', 'LineWidth', 2);
h4 = plot(x_plot_ppm, y_infu_plot, 'r-', 'LineWidth', 2);

xlabel('Offset (ppm)', 'FontSize', 10); 
ylabel('S/S_0', 'FontSize', 10);
if success_base && success_infu
    title(sprintf('B0-corrected (FWHM: %.2f/%.2f ppm)', FWHM_base, FWHM_infu), 'FontSize', 10);
else
    title('B0-corrected Z-spectra', 'FontSize', 10);
end
% leg = legend([h1 h3 h2 h4], {'Base (raw)', 'Base (fit)', 'Infu (raw)', 'Infu (fit)'}, ...
%     'Location', 'northeast', 'FontSize', 8);
% leg.ItemTokenSize = [15, 9];  % Reduce legend icon size
grid on;
set(gca, 'XDir', 'reverse', 'FontSize', 9); 
xlim([-4 4]); ylim([0 1.1]);

%% Subplot 2: After FWHM spillover correction
subplot(1,4,2);
y_infu_corrected_plot = ppval(pp_infu_corrected, x_plot_Hz);

h1 = plot(x_plot_ppm, y_base_plot, 'b-', 'LineWidth', 2); hold on;
h2 = plot(x_plot_ppm, y_infu_corrected_plot, 'Color', [1 0.5 0], 'LineWidth', 2);

if length(SelectPpm) == 1
    plot(SelectPpm, S_base_at_freq, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    plot(SelectPpm, S_infu_at_freq, 'o', 'Color', [1 0.5 0], 'MarkerSize', 10, 'LineWidth', 2);
end

xlabel('Offset (ppm)', 'FontSize', 10); 
ylabel('S/S_0', 'FontSize', 10);
if correction_applied && ~isnan(FWHM_base)
    title(sprintf('FWHM corrected (%.2f)', spillover_ratio), 'FontSize', 10);
elseif correction_applied
    title(sprintf('Depth-ratio corr (%.2f)', spillover_ratio), 'FontSize', 10);
else
    title('Uncorrected', 'FontSize', 10);
end
% leg = legend([h1 h2], {'Base', 'Infu'}, 'Location', 'northeast', 'FontSize', 8);
% leg.ItemTokenSize = [15, 9];
grid on;
xlim([-4 4]); ylim([0 1.1]); 
set(gca, 'XDir', 'reverse', 'FontSize', 9);

%% Subplot 3: Direct difference (Post - Pre)
subplot(1,4,3);
Delta_uncorrected = (y_infu_plot - y_base_plot) * 100;  % Post - Pre
Delta_corrected = (y_infu_corrected_plot - y_base_plot) * 100;  % Post - Pre

h1 = plot(x_plot_ppm, Delta_uncorrected, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); hold on;
h2 = plot(x_plot_ppm, smooth(Delta_corrected, 50), 'k-', 'LineWidth', 2);

if length(SelectPpm) == 1
    DeltaS_display = -DeltaS_Result.DeltaS_percent;  % Flip to Post - Pre
    plot(SelectPpm, DeltaS_display, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
    text(SelectPpm, DeltaS_display - 1, ...
        sprintf('  %.2f%%', DeltaS_display), 'FontSize', 9);
end

xlabel('Offset (ppm)', 'FontSize', 10); 
ylabel('\DeltaS (%)', 'FontSize', 10);
title('Direct diff (Post-Pre)', 'FontSize', 10);
% leg = legend([h1 h2], {'Raw', 'FWHM-corr'}, 'Location', 'best', 'FontSize', 8);
% leg.ItemTokenSize = [15, 9];
grid on;
xlim([-4 4]); 
set(gca, 'XDir', 'reverse', 'FontSize', 9);
xl = xlim;
plot(xl, [0 0], 'k:', 'LineWidth', 1);

%% Subplot 4: MTRasym analysis (Post - Pre)
subplot(1,4,4);

MTRasym_base = zeros(1, length(x_plot_ppm));
MTRasym_infu = zeros(1, length(x_plot_ppm));

for i = 1:length(x_plot_ppm)
    if x_plot_ppm(i) >= 0
        freq_pos = x_plot_ppm(i) * 400;
        freq_neg = -x_plot_ppm(i) * 400;
        
        S_neg_base = ppval(pp_base, freq_neg);
        S_pos_base = ppval(pp_base, freq_pos);
        MTRasym_base(i) = (S_neg_base - S_pos_base) * 100;
        
        S_neg_infu = ppval(pp_infu_corrected, freq_neg);
        S_pos_infu = ppval(pp_infu_corrected, freq_pos);
        MTRasym_infu(i) = (S_neg_infu - S_pos_infu) * 100;
    else
        MTRasym_base(i) = NaN;
        MTRasym_infu(i) = NaN;
    end
end

Delta_MTRasym = MTRasym_infu - MTRasym_base;  % Post - Pre

pos_idx = x_plot_ppm >= 0;
x_pos = x_plot_ppm(pos_idx);
MTRasym_base_pos = MTRasym_base(pos_idx);
MTRasym_infu_pos = MTRasym_infu(pos_idx);
Delta_MTRasym_pos = Delta_MTRasym(pos_idx);

h1 = plot(x_pos, MTRasym_base_pos, 'b-', 'LineWidth', 1.5); hold on;
h2 = plot(x_pos, MTRasym_infu_pos, 'Color', [1 0.5 0], 'LineWidth', 1.5);
h3 = plot(x_pos, smooth(Delta_MTRasym_pos, 30), 'k--', 'LineWidth', 2);

if length(SelectPpm) == 1 && SelectPpm > 0
    idx_sel = find(abs(x_pos - SelectPpm) < 0.05, 1);
    if ~isempty(idx_sel)
        plot(SelectPpm, Delta_MTRasym_pos(idx_sel), 'mo', 'MarkerSize', 10, 'LineWidth', 2);
        text(SelectPpm, Delta_MTRasym_pos(idx_sel) - 0.5, ...
            sprintf('  %.2f%%', Delta_MTRasym_pos(idx_sel)), 'FontSize', 9);
    end
end

xlabel('Offset (ppm)', 'FontSize', 10); 
ylabel('MTR_{asym} (%)', 'FontSize', 10);
title('Asymmetry (Post-Pre)', 'FontSize', 10);
% leg = legend([h1 h2 h3], {'Base', 'Infu', '\DeltaMTR_{asym}'}, 'Location', 'best', 'FontSize', 8);
% leg.ItemTokenSize = [15, 9];
grid on;
xlim([0 4]);
set(gca, 'FontSize', 9);
xl = xlim;
plot(xl, [0 0], 'k:', 'LineWidth', 1);

%% Figure title
if length(SelectPpm) == 1
    title_str = sprintf('%s - %s - %.2f ppm', AnimalID, ROI_name, SelectPpm);
else
    title_str = sprintf('%s - %s - %.1f-%.1f ppm', AnimalID, ROI_name, SelectPpm(1), SelectPpm(2));
end

annotation('textbox', [0 0.96 1 0.04], ...
    'String', title_str, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none');


%% Save results
DeltaS_Result.MTRasym_base = MTRasym_base;
DeltaS_Result.MTRasym_infu = MTRasym_infu;
DeltaS_Result.Delta_MTRasym = Delta_MTRasym;
DeltaS_Result.x_plot_ppm = x_plot_ppm;

% Export CSV (B0-corrected, no FWHM correction)
x_ppm_csv  = linspace(-4, 4, 4001)';
x_Hz_csv   = x_ppm_csv * 400;
y_base_csv = ppval(pp_base, x_Hz_csv);
y_infu_csv = ppval(pp_infu, x_Hz_csv);

T = table(x_ppm_csv, y_base_csv, y_infu_csv, ...
    'VariableNames', {'x_ppm', 'y_base', 'y_infu'});
csv_name = sprintf('Zspectrum_%s_%s.csv', AnimalID, ROI_name);
writetable(T, fullfile(pathBase, csv_name));
fprintf('CSV saved: %s\n', csv_name);


cd(pathBase);

if length(SelectPpm) == 1
    save_name = sprintf('DGE_%s_%s_%.2fppm_FWHMcorr.mat', AnimalID, ROI_name, SelectPpm);
    fig_name = sprintf('DGE_%s_%s_%.2fppm_FWHMcorr.fig', AnimalID, ROI_name, SelectPpm);
    png_name = sprintf('DGE_%s_%s_%.2fppm_FWHMcorr.png', AnimalID, ROI_name, SelectPpm);
else
    save_name = sprintf('DGE_%s_%s_%.1f-%.1fppm_FWHMcorr.mat', AnimalID, ROI_name, SelectPpm(1), SelectPpm(2));
    fig_name = sprintf('DGE_%s_%s_%.1f-%.1fppm_FWHMcorr.fig', AnimalID, ROI_name, SelectPpm(1), SelectPpm(2));
    png_name = sprintf('DGE_%s_%s_%.1f-%.1fppm_FWHMcorr.png', AnimalID, ROI_name, SelectPpm(1), SelectPpm(2));
end


save(save_name, 'DeltaS_Result');
saveas(gcf, fig_name);
saveas(gcf, png_name);
fprintf('\nResults saved to: %s\n', save_name);

end

%% Helper: Calculate water FWHM with fallback
function [fwhm_ppm, min_val, success] = calculate_water_FWHM(pp, Hz_per_ppm)
    x_search = linspace(-2, 2, 2000) * Hz_per_ppm;
    y_search = ppval(pp, x_search);
    
    [min_val, min_idx] = min(y_search);
    
    if min_val > 0.3
        fwhm_ppm = NaN;
        success = false;
        return;
    end
    
    half_max = (1 + min_val) / 2;
    
    idx_left = find(y_search(1:min_idx) > half_max, 1, 'last');
    idx_right = find(y_search(min_idx:end) > half_max, 1, 'first') + min_idx - 1;
    
    if isempty(idx_left) || isempty(idx_right)
        three_quarter_max = (1 + min_val) * 0.75;
        idx_left = find(y_search(1:min_idx) > three_quarter_max, 1, 'last');
        idx_right = find(y_search(min_idx:end) > three_quarter_max, 1, 'first') + min_idx - 1;
        
        if isempty(idx_left) || isempty(idx_right)
            fwhm_ppm = NaN;
            success = false;
        else
            fwhm_75 = (x_search(idx_right) - x_search(idx_left)) / Hz_per_ppm;
            fwhm_ppm = fwhm_75 * 1.5;
            success = true;
        end
    else
        fwhm_ppm = (x_search(idx_right) - x_search(idx_left)) / Hz_per_ppm;
        success = true;
    end
end
