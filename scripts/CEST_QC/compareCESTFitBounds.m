function cmp = compareCESTFitBounds(sessionDirs, varargin)
% COMPARECESTFITBOUNDS  Old-bounds vs new-bounds comparison of the pixelwise
%   Lorentzian fit for one animal (typically {CEST_base, CEST_infu}), run
%   BEFORE committing to a full 20+ h batch refit.
%
%   cmp = compareCESTFitBounds({baseDir, infuDir})
%
%   Each session directory must contain
%       LorentzianFitResult_oldBounds.mat   old Lb/Ub/P0   ('OldFitFile')
%       LorentzianFitResult.mat             new Lb/Ub/P0   ('NewFitFile')
%   plus the ROI files listed in CESTQCConfig().RoiFile.
%
%   For every session x variant the SAME downstream chain is run
%   (computeCESTQualityMetrics -> plotCESTGlucoseROIResults) with the
%   thresholds from CESTQCConfig, in a temporary folder, so nothing in the
%   session folders is touched.
%
%   COMPARED (per ROI, per session, old vs new)
%     pixelwise : median/IQR of glucose amplitude (17), width (18), position (19)
%                 fraction of ROI pixels with width or position AT a bound
%                 median RMSE_map, median R2_map, badFraction, nClean
%                 ROI mean of gluc_integral_fullMap (raw and post-QC)
%     ROI-mean  : popt_roi_post(17:19) and gluc_integral_roiFit_post
%     paired    : infu - base of the post-QC ROI mean and of the ROI-fit integral,
%                 for old and for new bounds
%
%   OUTPUT  cmp struct (also saved), <OutDir>/FitBoundsComparison_summary.txt,
%           <OutDir>/<roi>_FitBoundsComparison.png,
%           <OutDir>/<session>_<variant>_<roi>_MeanROI_Zspectrum_*.png
%
%   NAME-VALUE PAIRS
%     'SessionLabels' (default {'base','infu'})   one per sessionDirs entry
%     'OldFitFile'    (default 'LorentzianFitResult_oldBounds.mat')
%     'NewFitFile'    (default 'LorentzianFitResult.mat')
%     'RoiFile'       (default CESTQCConfig().RoiFile)
%     'GlucIdx'       (default [17 18 19])  amplitude, width, position
%     'BoundTol'      (default 1e-3)        relative tolerance for "at bound"
%     'OutDir'        (default <parent of first session>/FitBoundsComparison_<stamp>)
%     'ClimGl'        (default [0 0.02])    colour scale of integral maps
%
%   MATLAB R2018a compatible.

def = CESTQCConfig();
p = inputParser;
addParameter(p, 'SessionLabels', {'base', 'infu'});
addParameter(p, 'OldFitFile', 'LorentzianFitResult_oldBounds.mat');
addParameter(p, 'NewFitFile', 'LorentzianFitResult.mat');
addParameter(p, 'RoiFile', def.RoiFile);
addParameter(p, 'GlucIdx', [17 18 19]);
addParameter(p, 'BoundTol', 1e-3);
addParameter(p, 'OutDir', '');
addParameter(p, 'ClimGl', [0 0.02]);
parse(p, varargin{:});
opt = p.Results;

if ischar(sessionDirs), sessionDirs = {sessionDirs}; end
nSes = numel(sessionDirs);
if numel(opt.SessionLabels) ~= nSes
    error('compareCESTFitBounds:labels', 'SessionLabels must have one entry per session directory.');
end
roiFiles = cellstr(opt.RoiFile);
nROI = numel(roiFiles);
variants = {'old', 'new'};
fitFiles = {opt.OldFitFile, opt.NewFitFile};
iA = opt.GlucIdx(1); iW = opt.GlucIdx(2); iPos = opt.GlucIdx(3);

stamp = datestr(now, 'yyyymmdd_HHMMSS');
if isempty(opt.OutDir)
    opt.OutDir = fullfile(fileparts(sessionDirs{1}), sprintf('FitBoundsComparison_%s', stamp));
end
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
fprintf('Output -> %s\n', opt.OutDir);

%% ---------------- run QC + ROI fit for every session x variant ----------------
R = struct();   % R.(session).(variant)
tmpRoot = fullfile(tempdir, sprintf('compareCESTFitBounds_%s', stamp));
mkdir(tmpRoot);

for s = 1:nSes
    ses = opt.SessionLabels{s};
    for v = 1:2
        var = variants{v};
        srcFit = fullfile(sessionDirs{s}, fitFiles{v});
        if ~isfile(srcFit)
            error('compareCESTFitBounds:missingFit', 'Not found: %s', srcFit);
        end
        tmpDir = fullfile(tmpRoot, [ses '_' var]);
        mkdir(tmpDir);
        copyfile(srcFit, fullfile(tmpDir, 'LorentzianFitResult.mat'));
        for k = 1:nROI
            copyfile(fullfile(sessionDirs{s}, roiFiles{k}), fullfile(tmpDir, roiFiles{k}));
        end
        fprintf('\n===== %s / %s bounds =====\n', ses, var);

        F = load(srcFit, 'popt', 'P', 'gluc_integral_fullMap', 'Segment');
        [meanRMSE, meanR2, badFraction, roi_mean_clean, roi_std_clean, roi_labels, badPixelMask, qcFlags] = ...
            computeCESTQualityMetrics(tmpDir, 'RoiFile', roiFiles, 'SaveTxt', false, 'SaveMat', false);
        [~, ~, ~, roiFit] = plotCESTGlucoseROIResults(tmpDir, 'RoiFile', roiFiles, ...
            'BadPixelMask', badPixelMask, 'QCFlags', qcFlags, 'PlotMaps', false, ...
            'SavePostQCFit', false, 'SavePrefix', [ses '_' var '_']);
        pngs = dir(fullfile(tmpDir, '*.png'));
        for q = 1:numel(pngs)
            copyfile(fullfile(tmpDir, pngs(q).name), fullfile(opt.OutDir, pngs(q).name));
        end

        Lb = F.P.FIT.lower_limit_fit;  Ub = F.P.FIT.upper_limit_fit;  P0 = F.P.FIT.start_fit;
        ampMap = squeeze(F.popt(:,:,1,iA));
        widMap = squeeze(F.popt(:,:,1,iW));
        posMap = squeeze(F.popt(:,:,1,iPos));
        RMSE_map = qcFlags.values.RMSE;
        R2_map   = qcFlags.values.R2;
        tissueMask = (F.Segment(:,:,1) == 1);

        out = struct();
        out.fitFile = srcFit;
        out.Lb = Lb; out.Ub = Ub; out.P0 = P0;
        out.meanRMSE = meanRMSE; out.meanR2 = meanR2; out.badFraction = badFraction;
        out.roi_labels = roi_labels;
        out.gluc_integral_fullMap = F.gluc_integral_fullMap;
        out.badPixelMask = badPixelMask;
        out.tissueMask = tissueMask;
        out.roiFit = roiFit;
        for k = 1:nROI
            S = load(fullfile(sessionDirs{s}, roiFiles{k}));
            ROI_logical = logical(S.ROI);
            ROI_clean = ROI_logical & ~badPixelMask;
            out.roi(k).label = roi_labels{k};
            out.roi(k).ROI_logical = ROI_logical;
            out.roi(k).ROI_clean = ROI_clean;
            out.roi(k).nROI = sum(ROI_logical(:));
            out.roi(k).nClean = sum(ROI_clean(:));
            out.roi(k).amp = ampMap(ROI_clean);
            out.roi(k).wid = widMap(ROI_clean);
            out.roi(k).pos = posMap(ROI_clean);
            out.roi(k).integral = F.gluc_integral_fullMap(ROI_clean);
            out.roi(k).RMSE = RMSE_map(ROI_clean);
            out.roi(k).R2 = R2_map(ROI_clean);
            out.roi(k).fracWidAtBound = mean(local_atBound(out.roi(k).wid, Lb(iW), Ub(iW), opt.BoundTol));
            out.roi(k).fracPosAtBound = mean(local_atBound(out.roi(k).pos, Lb(iPos), Ub(iPos), opt.BoundTol));
            out.roi(k).fracAmpAtLb    = mean(abs(out.roi(k).amp - Lb(iA)) <= opt.BoundTol * max(1, abs(Lb(iA))));
            out.roi(k).roi_mean_raw   = mean(F.gluc_integral_fullMap(ROI_logical), 'omitnan');
            out.roi(k).roi_mean_clean = roi_mean_clean(k);
            out.roi(k).roi_std_clean  = roi_std_clean(k);
            out.roi(k).popt_roi_post  = roiFit(k).popt_roi_post;
            if isempty(out.roi(k).popt_roi_post)          % post-QC ROI fit skipped (no clean px)
                out.roi(k).popt_roi_post = nan(1, max(opt.GlucIdx));
            end
            out.roi(k).gluc_integral_roiFit_post = roiFit(k).gluc_integral_roiFit_post;
            out.roi(k).gluc_integral_roiFit_pre  = roiFit(k).gluc_integral_roiFit_pre;
        end
        R.(ses).(var) = out;
    end
end
rmdir(tmpRoot, 's');

%% ---------------- summary text ----------------
txtPath = fullfile(opt.OutDir, 'FitBoundsComparison_summary.txt');
fid = fopen(txtPath, 'w');
fprintf(fid, 'compareCESTFitBounds  %s\n', datestr(now));
for s = 1:nSes, fprintf(fid, '%s: %s\n', opt.SessionLabels{s}, sessionDirs{s}); end
fprintf(fid, 'QC: %s\n\n', qcFlags.thresholdString);

ses1 = opt.SessionLabels{1};
for v = 1:2
    o = R.(ses1).(variants{v});
    fprintf(fid, '[%s bounds]  glucose params %d/%d/%d (amp/width/pos)\n', variants{v}, iA, iW, iPos);
    fprintf(fid, '  Lb = %s\n  Ub = %s\n  P0 = %s\n', ...
        mat2str(o.Lb([iA iW iPos]), 4), mat2str(o.Ub([iA iW iPos]), 4), mat2str(o.P0([iA iW iPos]), 4));
end
fprintf(fid, '\n');

hdr = sprintf('%-8s %-5s %-9s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %10s %10s\n', ...
    'session', 'bnd', 'ROI', 'nROI', 'nClean', 'badFrac', 'medRMSE', 'medR2', ...
    'medAmp', 'medWid', 'medPos', 'fWid@B', 'fPos@B', 'fAmp@Lb', 'ROImeanCl', 'ROIfitInt', 'ROIfitWid');
fprintf(fid, '%s', hdr); fprintf('%s', hdr);
for k = 1:nROI
    for s = 1:nSes
        for v = 1:2
            o = R.(opt.SessionLabels{s}).(variants{v});
            r = o.roi(k);
            line = sprintf('%-8s %-5s %-9s %6d %6d %8.3f %8.4f %8.3f %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %10.5f %10.3f\n', ...
                opt.SessionLabels{s}, variants{v}, r.label, r.nROI, r.nClean, o.badFraction, ...
                median(r.RMSE, 'omitnan'), median(r.R2, 'omitnan'), ...
                median(r.amp), median(r.wid), median(r.pos), ...
                r.fracWidAtBound, r.fracPosAtBound, r.fracAmpAtLb, ...
                r.roi_mean_clean, r.gluc_integral_roiFit_post, r.popt_roi_post(iW));
            fprintf(fid, '%s', line); fprintf('%s', line);
        end
    end
end

% paired differences (only if exactly two sessions)
if nSes == 2
    sA = opt.SessionLabels{1}; sB = opt.SessionLabels{2};
    fprintf(fid, '\n[paired: %s - %s]\n', sB, sA);
    fprintf(fid, '%-9s %-5s %14s %14s %14s\n', 'ROI', 'bnd', 'dROImeanClean', 'dROIfitInt', 'ratioROImeanCl');
    for k = 1:nROI
        for v = 1:2
            a = R.(sA).(variants{v}).roi(k);  b = R.(sB).(variants{v}).roi(k);
            dMean = b.roi_mean_clean - a.roi_mean_clean;
            dFit  = b.gluc_integral_roiFit_post - a.gluc_integral_roiFit_post;
            ratio = b.roi_mean_clean / a.roi_mean_clean;
            line = sprintf('%-9s %-5s %14.5f %14.5f %14.3f\n', a.label, variants{v}, dMean, dFit, ratio);
            fprintf(fid, '%s', line); fprintf('%s', line);
        end
    end
end

% pixelwise agreement old vs new (same pixels: intersection of clean masks)
fprintf(fid, '\n[pixelwise old vs new, pixels clean in both]\n');
fprintf(fid, '%-8s %-9s %6s %8s %10s %10s\n', 'session', 'ROI', 'n', 'corr', 'medOldInt', 'medNewInt');
for s = 1:nSes
    ses = opt.SessionLabels{s};
    for k = 1:nROI
        both = R.(ses).old.roi(k).ROI_clean & R.(ses).new.roi(k).ROI_clean;
        xo = R.(ses).old.gluc_integral_fullMap(both);
        xn = R.(ses).new.gluc_integral_fullMap(both);
        ok = ~isnan(xo) & ~isnan(xn);
        c = NaN;
        if sum(ok) > 2, cc = corrcoef(xo(ok), xn(ok)); c = cc(1,2); end
        line = sprintf('%-8s %-9s %6d %8.3f %10.5f %10.5f\n', ses, R.(ses).old.roi(k).label, sum(ok), c, median(xo(ok)), median(xn(ok)));
        fprintf(fid, '%s', line); fprintf('%s', line);
    end
end
fclose(fid);
fprintf('\nSummary -> %s\n', txtPath);

%% ---------------- figures, one per ROI ----------------
for k = 1:nROI
    lbl = R.(ses1).old.roi(k).label;
    figK = figure('Color', 'white', 'Position', [50 50 1600 900], 'Visible', 'off');

    % row 1: integral maps inside ROI, old/new x sessions
    nCol = max(4, 2 * nSes);
    ci = 0;
    for v = 1:2
        for s = 1:nSes
            ci = ci + 1;
            o = R.(opt.SessionLabels{s}).(variants{v});
            m = o.gluc_integral_fullMap;
            m(~o.roi(k).ROI_clean) = NaN;
            subplot(3, nCol, ci);
            imagesc(rot90(m, 3), opt.ClimGl); axis image off; colormap(gca, jet);
            title(sprintf('%s / %s bounds (n=%d)', opt.SessionLabels{s}, variants{v}, o.roi(k).nClean), 'FontSize', 10);
        end
    end
    cb = colorbar('Position', [0.92 0.70 0.012 0.2]); cb.Label.String = 'gluc integral';

    % row 2: histograms width / position / amplitude, old vs new (all sessions pooled)
    fields = {'wid', 'pos', 'amp'};
    names  = {sprintf('width (param %d)', iW), sprintf('position (param %d)', iPos), sprintf('amplitude (param %d)', iA)};
    for f = 1:3
        subplot(3, nCol, nCol + f);
        hold on;
        for v = 1:2
            vals = [];
            for s = 1:nSes, vals = [vals; R.(opt.SessionLabels{s}).(variants{v}).roi(k).(fields{f})(:)]; end %#ok<AGROW>
            histogram(vals, 40, 'Normalization', 'probability', 'FaceAlpha', 0.45, 'DisplayName', [variants{v} ' bounds']);
        end
        hold off; box on; grid on;
        title(names{f}, 'FontSize', 10); legend('Location', 'best', 'FontSize', 8);
    end
    % scatter old vs new integral, pixels clean in both, all sessions
    subplot(3, nCol, nCol + 4);
    hold on;
    for s = 1:nSes
        ses = opt.SessionLabels{s};
        both = R.(ses).old.roi(k).ROI_clean & R.(ses).new.roi(k).ROI_clean;
        plot(R.(ses).old.gluc_integral_fullMap(both), R.(ses).new.gluc_integral_fullMap(both), '.', 'DisplayName', ses);
    end
    lim = [0 max(opt.ClimGl)]; plot(lim, lim, 'k--', 'HandleVisibility', 'off');
    hold off; axis square; box on; grid on;
    xlabel('integral, old bounds'); ylabel('integral, new bounds');
    title('pixelwise integral, clean in both', 'FontSize', 10); legend('Location', 'best', 'FontSize', 8);

    % row 3: bar charts ROI mean clean / ROI-fit integral / RMSE
    subplot(3, nCol, 2*nCol + 1);
    local_bars(R, opt.SessionLabels, variants, k, 'roi_mean_clean', 'ROI mean (post-QC)');
    subplot(3, nCol, 2*nCol + 2);
    local_bars(R, opt.SessionLabels, variants, k, 'gluc_integral_roiFit_post', 'ROI-mean fit integral (post-QC)');
    subplot(3, nCol, 2*nCol + 3);
    local_barsMedian(R, opt.SessionLabels, variants, k, 'RMSE', 'median RMSE (clean px)');
    subplot(3, nCol, 2*nCol + 4);
    local_barsFrac(R, opt.SessionLabels, variants, k);

    annotation(figK, 'textbox', [0 0.965 1 0.035], 'String', ...
        sprintf('%s  |  old vs new glucose-pool bounds  |  QC: %s', lbl, qcFlags.thresholdString), ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
    set(figK, 'PaperPositionMode', 'auto');
    print(figK, fullfile(opt.OutDir, sprintf('%s_FitBoundsComparison.png', lbl)), '-dpng', '-r200');
    close(figK);
end

cmp = struct('R', R, 'opt', opt, 'sessionDirs', {sessionDirs}, 'thresholdString', qcFlags.thresholdString);
save(fullfile(opt.OutDir, 'FitBoundsComparison.mat'), 'cmp', '-v7');
fprintf('Done -> %s\n', opt.OutDir);
end


%% ===================== local functions =====================

function tf = local_atBound(x, lb, ub, tol)
tf = (abs(x - lb) <= tol * max(1, abs(lb))) | (abs(x - ub) <= tol * max(1, abs(ub)));
end

function local_bars(R, sesLabels, variants, k, field, ttl)
nSes = numel(sesLabels);
vals = zeros(nSes, 2);
for s = 1:nSes
    for v = 1:2
        vals(s, v) = R.(sesLabels{s}).(variants{v}).roi(k).(field);
    end
end
bar(vals); set(gca, 'XTickLabel', sesLabels);
legend(variants, 'Location', 'best', 'FontSize', 8); title(ttl, 'FontSize', 10); grid on; box on;
end

function local_barsMedian(R, sesLabels, variants, k, field, ttl)
nSes = numel(sesLabels);
vals = zeros(nSes, 2);
for s = 1:nSes
    for v = 1:2
        vals(s, v) = median(R.(sesLabels{s}).(variants{v}).roi(k).(field), 'omitnan');
    end
end
bar(vals); set(gca, 'XTickLabel', sesLabels);
legend(variants, 'Location', 'best', 'FontSize', 8); title(ttl, 'FontSize', 10); grid on; box on;
end

function local_barsFrac(R, sesLabels, variants, k)
nSes = numel(sesLabels);
vals = zeros(nSes, 4);
names = cell(1, 4);
c = 0;
for v = 1:2
    c = c + 1; names{c} = [variants{v} ' width@bound'];
    for s = 1:nSes, vals(s, c) = R.(sesLabels{s}).(variants{v}).roi(k).fracWidAtBound; end
    c = c + 1; names{c} = [variants{v} ' pos@bound'];
    for s = 1:nSes, vals(s, c) = R.(sesLabels{s}).(variants{v}).roi(k).fracPosAtBound; end
end
bar(vals); set(gca, 'XTickLabel', sesLabels); ylim([0 1]);
legend(names, 'Location', 'best', 'FontSize', 8); title('fraction of clean px at a bound', 'FontSize', 10); grid on; box on;
end
