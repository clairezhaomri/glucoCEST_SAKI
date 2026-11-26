function [meanRMSE, meanR2, badFraction, roi_mean_clean, roi_std_clean, roi_labels, badPixelMask, qcFlags] = ...
    computeCESTQualityMetrics(dataDir, varargin)
% COMPUTECESTQUALITYMETRICS  Pixelwise fit-quality QC from an existing
%   LorentzianFitResult.mat, WITHOUT rerunning the fit.
%
%   [meanRMSE, meanR2, badFraction, roi_mean_clean, roi_std_clean, roi_labels, ...
%    badPixelMask, qcFlags] = computeCESTQualityMetrics(dataDir, ...)
%
%   ALL THRESHOLD DEFAULTS COME FROM CESTQCConfig.m. Edit thresholds there.
%   Name-value pairs passed here override the config for this call only and
%   are listed under "Overrides" in QC_summary.txt.
%
%   OUTPUTS (names identical to previous version)
%     meanRMSE        mean RMSE_map over tissueMask
%     meanR2          mean R2_map   over tissueMask
%     badFraction     sum(badPixelMask) / nTissue
%     roi_mean_clean  1 x nROI  mean gluc_integral_fullMap in ROI & ~badPixelMask
%     roi_std_clean   1 x nROI
%     roi_labels      1 x nROI cell, e.g. {'ctx','medulla'}
%     badPixelMask    logical [nx ny], union of all criteria
%     qcFlags         struct:
%       .masks.Zmin / .Zmax / .RMSE / .R2 / .Bound / .Spatial   logical [nx ny]
%       .values.Zmin / .Zmax / .RMSE / .R2                      double  [nx ny]
%       .thresholds.<every field of CESTQCConfig> + .RMSEThreshUsed
%       .thresholdString   one-line human-readable summary of thresholds
%       .overrides         cell of parameter names that differ from CESTQCConfig
%       .tissue            struct: nTissue, nFitted, per-criterion counts, nBad
%       .roi(k)            struct per ROI: label, nROI, nFitted, nBad, nClean,
%                          meanRaw, sdRaw, meanClean, sdClean, per-criterion counts
%
%   NAME-VALUE PAIRS (defaults: see CESTQCConfig.m)
%   'RoiFile', 'ZminThresh', 'ZmaxThresh', 'UseRMSEPercentile',
%   'RMSEPercentile', 'RMSEThresh', 'R2Thresh', 'UseBoundCheck', 'BoundTol',
%   'UseSpatialOutlier', 'SpatialMADMultiplier', 'SpatialWindow'
%   'SaveTxt'    (default true)   write <SavePrefix>QC_summary.txt into dataDir
%   'SaveMat'    (default true)   write <SavePrefix>QC_result.mat  into dataDir
%                                 (badPixelMask, qcFlags, opt; -v7, readable by
%                                 scipy.io.loadmat)
%   'SaveMaps'   (default false)  write <SavePrefix>badPixelMask.png
%   'SavePrefix' (default '')
%
%   MATLAB R2018a compatible.

def = CESTQCConfig();

p = inputParser;
addParameter(p, 'RoiFile',              def.RoiFile);
addParameter(p, 'ZminThresh',           def.ZminThresh);
addParameter(p, 'ZmaxThresh',           def.ZmaxThresh);
addParameter(p, 'UseRMSEPercentile',    def.UseRMSEPercentile);
addParameter(p, 'RMSEPercentile',       def.RMSEPercentile);
addParameter(p, 'RMSEThresh',           def.RMSEThresh);
addParameter(p, 'R2Thresh',             def.R2Thresh);
addParameter(p, 'UseBoundCheck',        def.UseBoundCheck);
addParameter(p, 'BoundTol',             def.BoundTol);
addParameter(p, 'UseSpatialOutlier',    def.UseSpatialOutlier);
addParameter(p, 'SpatialMADMultiplier', def.SpatialMADMultiplier);
addParameter(p, 'SpatialWindow',        def.SpatialWindow);
addParameter(p, 'SaveTxt',    true);
addParameter(p, 'SaveMat',    true);
addParameter(p, 'SaveMaps',   false);
addParameter(p, 'SavePrefix', '');
parse(p, varargin{:});
opt = p.Results;

% Which threshold parameters were explicitly overridden vs CESTQCConfig?
configNames = fieldnames(def);
overrides = setdiff(p.Parameters, p.UsingDefaults);
overrides = intersect(overrides, configNames);
overrides = overrides(:)';   % row cell
keepOv = true(size(overrides));
for q = 1:numel(overrides)
    if isequal(opt.(overrides{q}), def.(overrides{q})), keepOv(q) = false; end
end
overrides = overrides(keepOv);

%% -------- load fit result --------
fitResultsPath = fullfile(dataDir, 'LorentzianFitResult.mat');
if ~isfile(fitResultsPath)
    error('computeCESTQualityMetrics:missingFit', ...
        'LorentzianFitResult.mat not found in %s.', dataDir);
end
fitInfo = dir(fitResultsPath);
F = load(fitResultsPath, 'popt', 'P', 'gluc_integral_fullMap', 'M0_stack', 'Segment', 'Z_corrExt');
popt = F.popt; P = F.P; gluc_integral_fullMap = F.gluc_integral_fullMap;
Segment = F.Segment; Z_corrExt = F.Z_corrExt;

[nx, ny, ~, ~] = size(popt);
W = P.SEQ.w;
lb = P.FIT.lower_limit_fit;
ub = P.FIT.upper_limit_fit;

N_Zcorr = size(Z_corrExt, 4);
if numel(W) ~= N_Zcorr
    error('computeCESTQualityMetrics:offsetMismatch', ...
        ['P.SEQ.w has %d points but Z_corrExt has %d offsets in %s. ' ...
         'This LorentzianFitResult.mat is stale (built before the offset-count ' ...
         'fix). Delete it and re-run runCESTGlucoseROIFit for this subject.'], ...
        numel(W), N_Zcorr, dataDir);
end

%% -------- per-pixel maps --------
Zmin_map = squeeze(min(Z_corrExt(:,:,1,:), [], 4));
Zmax_map = squeeze(max(Z_corrExt(:,:,1,:), [], 4));

RMSE_map    = nan(nx, ny);
R2_map      = nan(nx, ny);
atBound_map = false(nx, ny);

tissueMask = (Segment(:,:,1) == 1);
tissueIdx  = find(tissueMask);
for idx = 1:numel(tissueIdx)
    [xi, yi] = ind2sub([nx, ny], tissueIdx(idx));

    measured = Z_corrExt(xi, yi, 1, :);
    measured = reshape(measured, [], 1);
    if any(isnan(measured))
        continue;
    end

    p_pixel = squeeze(popt(xi, yi, 1, :))';
    [f, ~] = fitmodelfunc_ANA(p_pixel', P);
    fitted = f(W);
    fitted = reshape(fitted, [], 1);

    if numel(measured) ~= numel(fitted)
        error('computeCESTQualityMetrics:pixelLengthMismatch', ...
            ['Pixel (%d,%d) in %s: measured has %d points but fitted ' ...
             'model output has %d points.'], ...
            xi, yi, dataDir, numel(measured), numel(fitted));
    end

    resid = measured - fitted;
    RMSE_map(xi, yi) = sqrt(mean(resid.^2));

    SStot = sum((measured - mean(measured)).^2);
    if SStot > 0
        R2_map(xi, yi) = 1 - sum(resid.^2) / SStot;
    end

    if any(abs(p_pixel - lb) < opt.BoundTol) || any(abs(p_pixel - ub) < opt.BoundTol)
        atBound_map(xi, yi) = true;
    end
end

nTissue = sum(tissueMask(:));
nFitted = sum(tissueMask(:) & ~isnan(RMSE_map(:)));

%% -------- criteria --------
if opt.UseRMSEPercentile
    rmseVals = RMSE_map(tissueMask);
    rmseThresh = prctile(rmseVals(~isnan(rmseVals)), opt.RMSEPercentile);
else
    rmseThresh = opt.RMSEThresh;
end
flagRMSE = (RMSE_map > rmseThresh) & tissueMask;

flagR2 = (R2_map < opt.R2Thresh) & tissueMask;

flagZmin = (Zmin_map > opt.ZminThresh) & tissueMask;
flagZmax = (Zmax_map > opt.ZmaxThresh) & tissueMask;

if opt.UseBoundCheck
    flagBound = atBound_map & tissueMask;
else
    flagBound = false(nx, ny);
end

if opt.UseSpatialOutlier
    localMed = medfilt2(gluc_integral_fullMap, opt.SpatialWindow);
    localDev = abs(gluc_integral_fullMap - localMed);
    madMap   = medfilt2(localDev, opt.SpatialWindow);
    flagSpatial = (localDev > opt.SpatialMADMultiplier * madMap) & tissueMask;
else
    flagSpatial = false(nx, ny);
end

badPixelMask = flagZmin | flagZmax | flagRMSE | flagR2 | flagBound | flagSpatial;

meanRMSE    = mean(RMSE_map(tissueMask), 'omitnan');
meanR2      = mean(R2_map(tissueMask), 'omitnan');
badFraction = sum(badPixelMask(:)) / max(nTissue, 1);

%% -------- qcFlags struct --------
qcFlags.masks.Zmin    = flagZmin;
qcFlags.masks.Zmax    = flagZmax;
qcFlags.masks.RMSE    = flagRMSE;
qcFlags.masks.R2      = flagR2;
qcFlags.masks.Bound   = flagBound;
qcFlags.masks.Spatial = flagSpatial;
qcFlags.values.Zmin   = Zmin_map;
qcFlags.values.Zmax   = Zmax_map;
qcFlags.values.RMSE   = RMSE_map;
qcFlags.values.R2     = R2_map;

for c = 1:numel(configNames)
    qcFlags.thresholds.(configNames{c}) = opt.(configNames{c});
end
qcFlags.thresholds.RMSEThreshUsed = rmseThresh;
qcFlags.overrides = overrides;

critNames = {'Zmin', 'Zmax', 'RMSE', 'R2', 'Bound', 'Spatial'};
qcFlags.tissue.nTissue = nTissue;
qcFlags.tissue.nFitted = nFitted;
qcFlags.tissue.nBad    = sum(badPixelMask(:));
for c = 1:numel(critNames)
    qcFlags.tissue.(critNames{c}) = sum(sum(qcFlags.masks.(critNames{c})));
end

qcFlags.thresholdString = local_thresholdString(opt, rmseThresh);

fprintf('  [QC criteria, tissue=%d px] %s\n', nTissue, qcFlags.thresholdString);
fprintf(['  [QC flags,    tissue=%d px] Zmin: %d (%.1f%%) | Zmax: %d (%.1f%%) | ' ...
         'RMSE: %d (%.1f%%) | R2: %d (%.1f%%) | Bound: %d (%.1f%%) | Spatial: %d (%.1f%%) | ' ...
         'bad(union): %d (%.1f%%)\n'], nTissue, ...
    qcFlags.tissue.Zmin,    100*qcFlags.tissue.Zmin/max(nTissue,1), ...
    qcFlags.tissue.Zmax,    100*qcFlags.tissue.Zmax/max(nTissue,1), ...
    qcFlags.tissue.RMSE,    100*qcFlags.tissue.RMSE/max(nTissue,1), ...
    qcFlags.tissue.R2,      100*qcFlags.tissue.R2/max(nTissue,1), ...
    qcFlags.tissue.Bound,   100*qcFlags.tissue.Bound/max(nTissue,1), ...
    qcFlags.tissue.Spatial, 100*qcFlags.tissue.Spatial/max(nTissue,1), ...
    qcFlags.tissue.nBad,    100*qcFlags.tissue.nBad/max(nTissue,1));

%% -------- per-ROI raw + cleaned stats --------
if ischar(opt.RoiFile) || isstring(opt.RoiFile)
    roiFiles = {char(opt.RoiFile)};
else
    roiFiles = cellstr(opt.RoiFile);
end

nROI = numel(roiFiles);
roi_mean_clean = nan(1, nROI);
roi_std_clean  = nan(1, nROI);
roi_labels     = cell(1, nROI);

roiInfo = struct('label', {}, 'file', {}, 'nROI', {}, 'nFitted', {}, 'nBad', {}, 'nClean', {}, ...
    'meanRaw', {}, 'sdRaw', {}, 'meanClean', {}, 'sdClean', {}, ...
    'Zmin', {}, 'Zmax', {}, 'RMSE', {}, 'R2', {}, 'Bound', {}, 'Spatial', {});

for k = 1:nROI
    roiFile = roiFiles{k};
    [~, roiLabel] = fileparts(roiFile);
    roiLabel = erase(roiLabel, 'ROI_');
    roi_labels{k} = roiLabel;

    roiInfo(k).label   = roiLabel;
    roiInfo(k).file    = roiFile;
    roiInfo(k).nROI    = NaN; roiInfo(k).nFitted = NaN;
    roiInfo(k).nBad    = NaN; roiInfo(k).nClean  = NaN;
    roiInfo(k).meanRaw = NaN; roiInfo(k).sdRaw   = NaN;
    roiInfo(k).meanClean = NaN; roiInfo(k).sdClean = NaN;
    for c = 1:numel(critNames)
        roiInfo(k).(critNames{c}) = NaN;
    end

    roiPath = fullfile(dataDir, roiFile);
    if ~isfile(roiPath)
        warning('computeCESTQualityMetrics:missingROI', 'ROI file not found, skipping: %s', roiPath);
        continue;
    end
    S = load(roiPath);
    if ~isfield(S, 'ROI')
        warning('computeCESTQualityMetrics:noROIvar', '%s has no variable ROI, skipping.', roiPath);
        continue;
    end
    ROI = S.ROI;

    if ~isequal(size(ROI), size(badPixelMask))
        error('computeCESTQualityMetrics:roiSizeMismatch', ...
            '%s is %s but badPixelMask is %s in %s.', ...
            roiFile, mat2str(size(ROI)), mat2str(size(badPixelMask)), dataDir);
    end

    ROI_logical = logical(ROI);
    ROI_clean   = ROI_logical & ~badPixelMask;

    croppedMap_raw = gluc_integral_fullMap;
    croppedMap_raw(~ROI_logical) = NaN;
    croppedMap_clean = gluc_integral_fullMap;
    croppedMap_clean(~ROI_clean) = NaN;

    roi_mean_clean(k) = mean(croppedMap_clean(:), 'omitnan');
    roi_std_clean(k)  = std(croppedMap_clean(:), 'omitnan');

    roiInfo(k).nROI      = sum(ROI_logical(:));
    roiInfo(k).nFitted   = sum(ROI_logical(:) & ~isnan(gluc_integral_fullMap(:)));
    roiInfo(k).nBad      = sum(ROI_logical(:) & badPixelMask(:));
    roiInfo(k).nClean    = sum(ROI_clean(:));
    roiInfo(k).meanRaw   = mean(croppedMap_raw(:), 'omitnan');
    roiInfo(k).sdRaw     = std(croppedMap_raw(:), 'omitnan');
    roiInfo(k).meanClean = roi_mean_clean(k);
    roiInfo(k).sdClean   = roi_std_clean(k);
    for c = 1:numel(critNames)
        roiInfo(k).(critNames{c}) = sum(sum(ROI_logical & qcFlags.masks.(critNames{c})));
    end
end
qcFlags.roi = roiInfo;

%% -------- save text summary (thresholds included) --------
if opt.SaveTxt
    txtPath = fullfile(dataDir, [opt.SavePrefix 'QC_summary.txt']);
    fid = fopen(txtPath, 'w');
    fprintf(fid, 'CEST fit quality summary\n');
    fprintf(fid, 'Folder: %s\n', dataDir);
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'FitResultFile: %s (modified %s)\n', fitResultsPath, datestr(fitInfo.datenum));
    fprintf(fid, 'QCFunction: computeCESTQualityMetrics | ConfigFile: %s\n\n', which('CESTQCConfig'));

    % --- legacy lines kept EXACTLY as before (safe for existing parsers) ---
    fprintf(fid, 'Mean RMSE (tissue): %.6f\n', meanRMSE);
    fprintf(fid, 'Mean R2 (tissue): %.6f\n', meanR2);
    fprintf(fid, 'Bad pixel fraction (tissue): %.4f\n\n', badFraction);
    for k = 1:nROI
        fprintf(fid, 'ROI: %s | Mean: %.4f | SD: %.4f\n', roi_labels{k}, roi_mean_clean(k), roi_std_clean(k));
    end

    % --- QC criteria actually used ---
    fprintf(fid, '\n[QC criteria]\n');
    fprintf(fid, 'ZminThresh: %g            (flag if Zmin_map > ZminThresh)\n', opt.ZminThresh);
    fprintf(fid, 'ZmaxThresh: %g            (flag if Zmax_map > ZmaxThresh)\n', opt.ZmaxThresh);
    fprintf(fid, 'UseRMSEPercentile: %d\n', opt.UseRMSEPercentile);
    fprintf(fid, 'RMSEPercentile: %g\n', opt.RMSEPercentile);
    fprintf(fid, 'RMSEThresh: %g\n', opt.RMSEThresh);
    fprintf(fid, 'RMSEThreshUsed: %.6g       (flag if RMSE_map > RMSEThreshUsed)\n', rmseThresh);
    fprintf(fid, 'R2Thresh: %g               (flag if R2_map < R2Thresh)\n', opt.R2Thresh);
    fprintf(fid, 'UseBoundCheck: %d\n', opt.UseBoundCheck);
    fprintf(fid, 'BoundTol: %g\n', opt.BoundTol);
    fprintf(fid, 'UseSpatialOutlier: %d\n', opt.UseSpatialOutlier);
    fprintf(fid, 'SpatialMADMultiplier: %g\n', opt.SpatialMADMultiplier);
    fprintf(fid, 'SpatialWindow: %s\n', mat2str(opt.SpatialWindow));
    fprintf(fid, 'RoiFile: %s\n', strjoin(roiFiles, ', '));
    if isempty(overrides)
        fprintf(fid, 'Overrides vs CESTQCConfig: none\n');
    else
        fprintf(fid, 'Overrides vs CESTQCConfig: %s\n', strjoin(overrides, ', '));
    end

    % --- tissue-level flag counts ---
    fprintf(fid, '\n[Tissue flags]\n');
    fprintf(fid, 'nTissue: %d\n', nTissue);
    fprintf(fid, 'nFitted: %d\n', nFitted);
    for c = 1:numel(critNames)
        fprintf(fid, '%s: %d (%.1f%%)\n', critNames{c}, qcFlags.tissue.(critNames{c}), ...
            100*qcFlags.tissue.(critNames{c})/max(nTissue,1));
    end
    fprintf(fid, 'nBad (union): %d (%.1f%%)\n', qcFlags.tissue.nBad, 100*qcFlags.tissue.nBad/max(nTissue,1));

    % --- per-ROI detail ---
    fprintf(fid, '\n[ROI detail]\n');
    for k = 1:nROI
        r = roiInfo(k);
        fprintf(fid, 'ROI: %s | nROI: %d | nFitted: %d | nBad: %d | nClean: %d | MeanRaw: %.4f | SDRaw: %.4f | MeanClean: %.4f | SDClean: %.4f\n', ...
            r.label, r.nROI, r.nFitted, r.nBad, r.nClean, r.meanRaw, r.sdRaw, r.meanClean, r.sdClean);
        fprintf(fid, 'ROI: %s flags | Zmin: %d | Zmax: %d | RMSE: %d | R2: %d | Bound: %d | Spatial: %d\n', ...
            r.label, r.Zmin, r.Zmax, r.RMSE, r.R2, r.Bound, r.Spatial);
    end
    fclose(fid);
    fprintf('QC summary saved: %s\n', txtPath);
end

%% -------- save masks + thresholds as .mat (retrospective / Python) --------
if opt.SaveMat
    matPath = fullfile(dataDir, [opt.SavePrefix 'QC_result.mat']);
    save(matPath, 'badPixelMask', 'qcFlags', 'opt', '-v7');
    fprintf('QC result saved: %s\n', matPath);
end

%% -------- optional bad-pixel overlay map --------
if opt.SaveMaps
    figQC = figure('Position', [100, 100, 500, 500], 'Visible', 'off');
    imagesc(F.M0_stack, [0, 150]);
    colormap gray
    axis image
    axis off
    hold on
    badRed = cat(3, ones(nx, ny), zeros(nx, ny), zeros(nx, ny));
    hBad = imshow(badRed);
    set(hBad, 'AlphaData', 0.5 * double(badPixelMask));
    title('Excluded pixels', 'FontSize', 16);
    hold off
    set(figQC, 'PaperPositionMode', 'auto');
    print(figQC, fullfile(dataDir, [opt.SavePrefix 'badPixelMask.png']), '-dpng', '-r300');
    close(figQC);
end

end


function s = local_thresholdString(opt, rmseThresh)
% One-line, human-readable summary of the thresholds actually applied.
if opt.UseRMSEPercentile
    rmsePart = sprintf('RMSE>%.4g (p%g)', rmseThresh, opt.RMSEPercentile);
else
    rmsePart = sprintf('RMSE>%.4g (abs)', rmseThresh);
end
if opt.UseBoundCheck
    boundPart = sprintf('Bound tol %.1g', opt.BoundTol);
else
    boundPart = 'Bound off';
end
if opt.UseSpatialOutlier
    spatialPart = sprintf('Spatial x%g [%dx%d]', opt.SpatialMADMultiplier, ...
        opt.SpatialWindow(1), opt.SpatialWindow(2));
else
    spatialPart = 'Spatial off';
end
s = sprintf('Zmin>%g | Zmax>%g | %s | R2<%g | %s | %s', ...
    opt.ZminThresh, opt.ZmaxThresh, rmsePart, opt.R2Thresh, boundPart, spatialPart);
end
