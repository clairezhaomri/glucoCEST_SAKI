function [roi_mean_before, roi_std_before, roi_mean_after, roi_std_after, roi_labels] = ...
    plotCESTQCBeforeAfter(dataDir, varargin)
% PLOTCESTQCBEFOREAFTER  Before/after QC comparison: glucose-overlay maps
%   PLUS example Z-spectra underneath each map.
%
%   [roi_mean_before, roi_std_before, roi_mean_after, roi_std_after, roi_labels] = ...
%       plotCESTQCBeforeAfter(dataDir, ...)
%
%   THRESHOLDS: this function does NOT define any QC threshold. It calls
%   computeCESTQualityMetrics, whose defaults come from CESTQCConfig.m.
%   So interactive plotting and WZ_QueueCEST_QC.m use identical thresholds
%   unless you deliberately pass 'QCOptions' here (interactive experiments).
%
%   TWO WAYS TO CALL
%   (a) interactive, recompute QC:
%         plotCESTQCBeforeAfter(dataDir)
%         plotCESTQCBeforeAfter(dataDir, 'QCOptions', {'ZminThresh', 0.2})
%   (b) batch, reuse masks already computed by computeCESTQualityMetrics:
%         plotCESTQCBeforeAfter(dataDir, 'BadPixelMask', badPixelMask, 'QCFlags', qcFlags)
%       Pass BOTH; without QCFlags the per-criterion reasons are unavailable.
%
%   NAME-VALUE PAIRS
%   'RoiFile'         (default CESTQCConfig().RoiFile) - string OR cell array
%   'RotN'            (default 3)
%   'ClimBg'          (default [0, 150])
%   'ClimGl'          (default [0, 0.2])
%   'AlphaGlucose'    (default 0.7)
%   'FigPosition'     (default [100 100 1400 750])
%   'SavePrefix'      (default '')
%   'NExampleSpectra' (default 3)
%   'FontSize'        (default 14)
%   'ShowThresholds'  (default true)  print qcFlags.thresholdString on the figure
%   'BadPixelMask'    (default [])
%   'QCFlags'         (default [])
%   'QCOptions'       (default {})   forwarded to computeCESTQualityMetrics
%
%   MATLAB R2018a compatible (no sgtitle / tiledlayout / exportgraphics).

def = CESTQCConfig();

p = inputParser;
addParameter(p, 'RoiFile', def.RoiFile);
addParameter(p, 'RotN', 3);
addParameter(p, 'ClimBg', [0, 150]);
addParameter(p, 'ClimGl', [0, 0.2]);
addParameter(p, 'AlphaGlucose', 0.7);
addParameter(p, 'FigPosition', [100, 100, 1400, 750]);
addParameter(p, 'SavePrefix', '');
addParameter(p, 'NExampleSpectra', 3);
addParameter(p, 'FontSize', 14);
addParameter(p, 'ShowThresholds', true);
addParameter(p, 'BadPixelMask', []);
addParameter(p, 'QCFlags', []);
addParameter(p, 'QCOptions', {});
parse(p, varargin{:});
opt = p.Results;

fitResultsPath = fullfile(dataDir, 'LorentzianFitResult.mat');
if ~isfile(fitResultsPath)
    error('plotCESTQCBeforeAfter:missingFit', ...
        'LorentzianFitResult.mat not found in %s. Run runCESTGlucoseROIFit first.', dataDir);
end
F = load(fitResultsPath, 'gluc_integral_fullMap', 'M0_stack', 'popt', 'P', 'Z_corrExt');
gluc_integral_fullMap = F.gluc_integral_fullMap;
M0_stack = F.M0_stack;
popt = F.popt;
P = F.P;
Z_corrExt = F.Z_corrExt;
W = P.SEQ.w;
[nx, ny, ~, ~] = size(popt);

if ischar(opt.RoiFile) || isstring(opt.RoiFile)
    roiFiles = {char(opt.RoiFile)};
else
    roiFiles = cellstr(opt.RoiFile);
end

if isempty(opt.BadPixelMask)
    % SaveTxt/SaveMat = false: interactive plotting must never overwrite the
    % QC_summary.txt / QC_result.mat written by the batch queue.
    [~, ~, ~, ~, ~, ~, badPixelMask, qcFlags] = computeCESTQualityMetrics(dataDir, ...
        'RoiFile', roiFiles, 'SaveTxt', false, 'SaveMat', false, opt.QCOptions{:});
else
    badPixelMask = opt.BadPixelMask;
    qcFlags = opt.QCFlags;
    if isempty(qcFlags)
        warning('plotCESTQCBeforeAfter:noQCFlags', ...
            ['BadPixelMask was given without QCFlags: per-criterion reasons and ' ...
             'threshold annotation are unavailable. Pass ''QCFlags'', qcFlags too.']);
    end
end

if ~isempty(qcFlags) && isfield(qcFlags, 'thresholdString')
    thresholdString = qcFlags.thresholdString;
else
    thresholdString = 'thresholds: n/a (qcFlags not provided)';
end

nROI = numel(roiFiles);
roi_mean_before = nan(1, nROI);
roi_std_before = nan(1, nROI);
roi_mean_after = nan(1, nROI);
roi_std_after = nan(1, nROI);
roi_labels = cell(1, nROI);

M0_rot = rot90(M0_stack, opt.RotN);
gluc_rot = rot90(gluc_integral_fullMap, opt.RotN);
badPixelMask_rot = rot90(badPixelMask, opt.RotN);

[rowGrid, colGrid] = ndgrid(1:nx, 1:ny);
rowGrid_rot = rot90(rowGrid, opt.RotN);
colGrid_rot = rot90(colGrid, opt.RotN);

nEx = opt.NExampleSpectra;
gridCols = 2 * nEx;
leftCols = 1:nEx;
rightCols = (nEx+1):(2*nEx);

for k = 1:nROI
    roiFile = roiFiles{k};
    [~, roiLabel] = fileparts(roiFile);
    roiLabel = erase(roiLabel, 'ROI_');
    roi_labels{k} = roiLabel;
    prefix = [opt.SavePrefix roiLabel '_'];

    roiPath = fullfile(dataDir, roiFile);
    if ~isfile(roiPath)
        warning('plotCESTQCBeforeAfter:missingROI', 'ROI file not found, skipping: %s', roiPath);
        continue;
    end
    S = load(roiPath);
    if ~isfield(S, 'ROI')
        warning('plotCESTQCBeforeAfter:noROIvar', '%s has no variable ROI, skipping.', roiPath);
        continue;
    end
    ROI = S.ROI;
    ROI_logical = logical(ROI);

    if ~isempty(qcFlags)
        nInROI = sum(sum(ROI_logical));
        nInROIandFitted = sum(sum(ROI_logical & ~isnan(gluc_integral_fullMap)));
        fprintf(['  [%s, this-ROI only, %d px total, %d actually fitted] Zmin: %d | Zmax: %d | RMSE: %d | R2: %d | ' ...
                 'Bound: %d | Spatial: %d\n'], ...
            roiLabel, nInROI, nInROIandFitted, ...
            sum(sum(ROI_logical & qcFlags.masks.Zmin)), ...
            sum(sum(ROI_logical & qcFlags.masks.Zmax)), ...
            sum(sum(ROI_logical & qcFlags.masks.RMSE)), ...
            sum(sum(ROI_logical & qcFlags.masks.R2)), ...
            sum(sum(ROI_logical & qcFlags.masks.Bound)), ...
            sum(sum(ROI_logical & qcFlags.masks.Spatial)));
    end

    ROI_rot = rot90(ROI, opt.RotN);
    ROI_nan_rot = double(ROI_rot);
    ROI_nan_rot(ROI_nan_rot == 0) = NaN;

    croppedMap_before = gluc_rot .* ROI_nan_rot;
    croppedMap_after = croppedMap_before;
    croppedMap_after(badPixelMask_rot) = NaN;

    roi_mean_before(k) = mean(croppedMap_before(:), 'omitnan');
    roi_std_before(k) = std(croppedMap_before(:), 'omitnan');
    roi_mean_after(k) = mean(croppedMap_after(:), 'omitnan');
    roi_std_after(k) = std(croppedMap_after(:), 'omitnan');

    fprintf('[%s] Before QC: Mean %.4f, SD %.4f | After QC: Mean %.4f, SD %.4f\n', ...
        roiLabel, roi_mean_before(k), roi_std_before(k), roi_mean_after(k), roi_std_after(k));

    goodPick = local_pickGoodExamples(ROI_logical, badPixelMask, qcFlags, gluc_integral_fullMap, nEx);
    [badPick, critNames] = local_pickOnePerCriterion(ROI_logical, badPixelMask, qcFlags, gluc_integral_fullMap, nEx);

    %% -------- figure: 2 maps (rows 1-2) + 2*nEx spectra (row 3) --------
    figBA = figure('Position', opt.FigPosition, 'Visible', 'off');

    row1idx = leftCols;
    row2idx = leftCols + gridCols;
    beforeCells = [row1idx, row2idx];

    row1idxR = rightCols;
    row2idxR = rightCols + gridCols;
    afterCells = [row1idxR, row2idxR];

    %% -------- Before QC panel --------
    ax_bg1 = subplot(3, gridCols, beforeCells);
    axPos1 = get(ax_bg1, 'Position');
    delete(ax_bg1);

    ax1 = axes('Position', axPos1);
    axes(ax1);
    imagesc(M0_rot, opt.ClimBg);
    colormap(ax1, gray);
    axis image;
    axis off;

    ax2 = axes('Position', axPos1);
    axes(ax2);
    h_overlay1 = imagesc(croppedMap_before);
    colormap(ax2, jet);
    axis image;
    set(h_overlay1, 'AlphaData', opt.AlphaGlucose * double(~isnan(croppedMap_before)));
    caxis(ax2, opt.ClimGl);
    ax2.XColor = 'none';
    ax2.YColor = 'none';
    ax2.Color = 'none';
    title(ax2, sprintf('Before QC - %s\nMean %.4f, SD %.4f', roiLabel, roi_mean_before(k), roi_std_before(k)), ...
        'FontSize', opt.FontSize);
    linkaxes([ax1, ax2]);
    ax1.Position = ax2.Position;

    %% -------- After QC panel (colorbar added BEFORE position sync) --------
    ax_bg2 = subplot(3, gridCols, afterCells);
    axPos2 = get(ax_bg2, 'Position');
    delete(ax_bg2);

    ax3 = axes('Position', axPos2);
    axes(ax3);
    imagesc(M0_rot, opt.ClimBg);
    colormap(ax3, gray);
    axis image;
    axis off;

    ax4 = axes('Position', axPos2);
    axes(ax4);
    h_overlay2 = imagesc(croppedMap_after);
    colormap(ax4, jet);
    axis image;
    set(h_overlay2, 'AlphaData', opt.AlphaGlucose * double(~isnan(croppedMap_after)));
    caxis(ax4, opt.ClimGl);
    ax4.XColor = 'none';
    ax4.YColor = 'none';
    ax4.Color = 'none';
    title(ax4, sprintf('After QC - %s\nMean %.4f, SD %.4f', roiLabel, roi_mean_after(k), roi_std_after(k)), ...
        'FontSize', opt.FontSize);

    cb = colorbar(ax4);
    cb.FontSize = opt.FontSize;
    cb.Label.String = 'Glucose integral';
    cb.Label.FontSize = opt.FontSize;

    linkaxes([ax3, ax4]);
    ax3.Position = ax4.Position;

    % --- Bottom row: nEx bad-pixel spectra (left) + nEx accepted (right) ---
    for e = 1:nEx
        cellIdxBad = 2*gridCols + e;
        subplot(3, gridCols, cellIdxBad);
        local_plotSpectrum(badPick(e), P, W, popt, Z_corrExt, rowGrid_rot, colGrid_rot, ...
            nx, ny, 'Bad', opt.FontSize, gluc_integral_fullMap, qcFlags, critNames{e});

        cellIdxGood = 2*gridCols + nEx + e;
        subplot(3, gridCols, cellIdxGood);
        local_plotSpectrum(goodPick(e), P, W, popt, Z_corrExt, rowGrid_rot, colGrid_rot, ...
            nx, ny, 'OK', opt.FontSize, gluc_integral_fullMap, qcFlags, '');
    end

    % --- thresholds actually used, printed on the figure (retrospective) ---
    if opt.ShowThresholds
        annotation(figBA, 'textbox', [0, 0.002, 1, 0.035], ...
            'String', ['QC: ' thresholdString], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'EdgeColor', 'none', 'Interpreter', 'none', ...
            'FontSize', max(8, round(opt.FontSize * 0.65)));
    end

    set(figBA, 'PaperPositionMode', 'auto');
    print(figBA, fullfile(dataDir, [prefix 'QC_before_after.png']), '-dpng', '-r300');
    close(figBA);
end

end


function picks = local_pickGoodExamples(ROI_logical, badPixelMask, qcFlags, gluc_integral_fullMap, nEx)
% Ranks accepted (non-bad) pixels by HIGHEST glucose integral first, then
% LOWEST RMSE, then HIGHEST R2. Explicitly excludes any pixel with NaN
% glucose/RMSE/R2 -- these are ROI pixels that fall outside the actual
% fitted tissue mask and were never evaluated by any QC criterion, so
% they must never be selected as "accepted" examples.
picks = nan(1, nEx);
goodMask = ROI_logical & ~badPixelMask & ~isnan(gluc_integral_fullMap);
candidateLinearIdx = find(goodMask);
n = numel(candidateLinearIdx);
if n == 0
    return;
end

gVals = gluc_integral_fullMap(candidateLinearIdx);

if ~isempty(qcFlags)
    rmseVals = qcFlags.values.RMSE(candidateLinearIdx);
    r2Vals = qcFlags.values.R2(candidateLinearIdx);
    validRow = ~isnan(rmseVals) & ~isnan(r2Vals);
    candidateLinearIdx = candidateLinearIdx(validRow);
    gVals = gVals(validRow);
    rmseVals = rmseVals(validRow);
    r2Vals = r2Vals(validRow);
    if isempty(candidateLinearIdx)
        return;
    end
    M = [gVals(:), rmseVals(:), r2Vals(:), candidateLinearIdx(:)];
    M = sortrows(M, [-1, 2, -3]);
    sortedIdx = M(:, 4);
else
    [~, sortIdx] = sort(gVals, 'descend');
    sortedIdx = candidateLinearIdx(sortIdx);
end

take = min(nEx, numel(sortedIdx));
picks(1:take) = sortedIdx(1:take)';
end


function [picks, critNames] = local_pickOnePerCriterion(ROI_logical, badPixelMask, qcFlags, gluc_integral_fullMap, nEx)
picks = nan(1, nEx);
critNames = repmat({'n/a'}, 1, nEx);

if isempty(qcFlags)
    badInROI = find(ROI_logical & badPixelMask & ~isnan(gluc_integral_fullMap));
    picks = local_pickWorst_(badInROI, gluc_integral_fullMap, nEx);
    critNames = repmat({'combined'}, 1, nEx);
    return;
end

priorityOrder = {'Zmin', 'Zmax', 'RMSE', 'R2', 'Bound', 'Spatial'};
activeCriteria = {};
for c = 1:numel(priorityOrder)
    name = priorityOrder{c};
    poolMask = ROI_logical & qcFlags.masks.(name) & ~isnan(gluc_integral_fullMap);
    if any(poolMask(:))
        activeCriteria{end+1} = name; %#ok<AGROW>
    end
end

if isempty(activeCriteria)
    return;
end

for e = 1:nEx
    name = activeCriteria{mod(e-1, numel(activeCriteria)) + 1};
    poolMask = ROI_logical & qcFlags.masks.(name) & ~isnan(gluc_integral_fullMap);
    candidateLinearIdx = find(poolMask);
    if isempty(candidateLinearIdx)
        continue;
    end
    vals = gluc_integral_fullMap(candidateLinearIdx);
    [~, sortIdx] = sort(vals, 'descend');
    sortedIdx = candidateLinearIdx(sortIdx);

    chosen = NaN;
    for cand = 1:numel(sortedIdx)
        if ~any(picks(1:e-1) == sortedIdx(cand))
            chosen = sortedIdx(cand);
            break;
        end
    end
    if ~isnan(chosen)
        picks(e) = chosen;
        critNames{e} = name;
    end
end
end


function picks = local_pickWorst_(candidateLinearIdx, valueMap, nEx)
picks = nan(1, nEx);
n = numel(candidateLinearIdx);
if n == 0
    return;
end
vals = valueMap(candidateLinearIdx);
[~, sortIdx] = sort(vals, 'descend');
sortedIdx = candidateLinearIdx(sortIdx);
take = min(nEx, n);
picks(1:take) = sortedIdx(1:take)';
end


function reasonStr = local_reasonString(xi, yi, qcFlags)
if isempty(qcFlags)
    reasonStr = 'n/a';
    return;
end
parts = {};
if qcFlags.masks.Zmin(xi, yi)
    parts{end+1} = sprintf('Zmin=%.3f>%.3g', qcFlags.values.Zmin(xi,yi), qcFlags.thresholds.ZminThresh);
end
if qcFlags.masks.Zmax(xi, yi)
    parts{end+1} = sprintf('Zmax=%.3f>%.3g', qcFlags.values.Zmax(xi,yi), qcFlags.thresholds.ZmaxThresh);
end
if qcFlags.masks.RMSE(xi, yi)
    parts{end+1} = sprintf('RMSE=%.4f>%.4g', qcFlags.values.RMSE(xi,yi), qcFlags.thresholds.RMSEThreshUsed);
end
if qcFlags.masks.R2(xi, yi)
    parts{end+1} = sprintf('R2=%.3f<%.3g', qcFlags.values.R2(xi,yi), qcFlags.thresholds.R2Thresh);
end
if qcFlags.masks.Bound(xi, yi)
    parts{end+1} = 'atBound';
end
if qcFlags.masks.Spatial(xi, yi)
    parts{end+1} = 'spatialOutlier';
end
if isempty(parts)
    reasonStr = 'unknown';
else
    reasonStr = strjoin(parts, ', ');
end
end


function acceptStr = local_acceptString(xi, yi, qcFlags)
if isempty(qcFlags)
    acceptStr = 'n/a';
    return;
end
acceptStr = sprintf('RMSE=%.4f, R2=%.3f', qcFlags.values.RMSE(xi,yi), qcFlags.values.R2(xi,yi));
end


function local_plotSpectrum(linIdx, P, W, popt, Z_corrExt, rowGrid_rot, colGrid_rot, nx, ny, label, fontSize, gluc_integral_fullMap, qcFlags, targetCrit)
if isnan(linIdx)
    axis off;
    title(sprintf('%s: n/a', label), 'FontSize', fontSize);
    return;
end
[xi, yi] = ind2sub([nx, ny], linIdx);

measured = reshape(Z_corrExt(xi, yi, 1, :), [], 1);
p_pixel = squeeze(popt(xi, yi, 1, :))';
[f, ~] = fitmodelfunc_ANA(p_pixel', P);
w_plot = linspace(min(W), max(W), 200)';

plot(W, measured, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
hold on;
plot(w_plot, f(w_plot), '-r', 'LineWidth', 1.5);
hold off;
set(gca, 'XDir', 'reverse', 'FontSize', fontSize);
xlim([min(W) max(W)]);
ylim([0 1]);
grid on;

matchMask = (rowGrid_rot == xi) & (colGrid_rot == yi);
[rr, rc] = find(matchMask, 1);
glucVal = gluc_integral_fullMap(xi, yi);
if isempty(rr)
    coordStr = '(?,?)';
else
    coordStr = sprintf('(%d,%d)', rr, rc);
end

if strcmp(label, 'Bad')
    reasonStr = local_reasonString(xi, yi, qcFlags);
    fprintf('  [Bad example - target %s] pixel %s g=%.3f -- flagged by: %s\n', ...
        targetCrit, coordStr, glucVal, reasonStr);
    title(sprintf('%s[%s] %s g=%.3f\n%s', label, targetCrit, coordStr, glucVal, reasonStr), ...
        'FontSize', fontSize*0.7);
else
    acceptStr = local_acceptString(xi, yi, qcFlags);
    fprintf('  [OK example] pixel %s g=%.3f -- accepted, %s\n', coordStr, glucVal, acceptStr);
    title(sprintf('%s %s g=%.3f\n%s', label, coordStr, glucVal, acceptStr), 'FontSize', fontSize*0.8);
end
end
