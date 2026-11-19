function [roi_mean, roi_std, roi_labels, roiFit] = plotCESTGlucoseROIResults(dataDir, varargin)
% PLOTCESTGLUCOSEROIRESULTS  (Re)plot mask/glucose overlay/Z-spectrum figures
%   from a previously saved LorentzianFitResult.mat, WITHOUT rerunning the
%   pixelwise fit. QC-AWARE: if a badPixelMask is available, the mean ROI
%   Z-spectrum is ALSO computed over ROI & ~badPixelMask, fitted again, and
%   compared with the pre-QC fit. The post-QC ROI fit is saved to
%   dataDir/LorentzianFitPostQC.mat.
%
%   [roi_mean, roi_std, roi_labels, roiFit] = plotCESTGlucoseROIResults(dataDir, ...)
%
%   OUTPUTS
%     roi_mean, roi_std  1 x nROI  PRE-QC pixelwise ROI stats of
%                        gluc_integral_fullMap (unchanged meaning, so
%                        WZ_QueueCEST_Analyze.m results stay comparable)
%     roi_labels         1 x nROI cell, e.g. {'ctx','medulla'}
%     roiFit             1 x nROI struct with the mean-ROI-spectrum fits:
%        .label .file .W
%        .nPix_pre  .meanZspec_pre  .popt_roi_pre  .gluc_integral_roiFit_pre
%        .nPix_post .meanZspec_post .popt_roi_post .gluc_integral_roiFit_post
%        .roi_mean_pre .roi_std_pre .roi_mean_post .roi_std_post  (pixelwise map stats)
%        .qcApplied (logical)
%       gluc_integral_roiFit_* = popt_roi(17) * popt_roi(18) * pi/2, i.e. the
%       same definition as gluc_integral_fullMap but from the ROI-mean fit.
%
%   WHERE THE QC MASK COMES FROM (first match wins)
%     1. 'BadPixelMask' (+ 'QCFlags') passed explicitly  (WZ_QueueCEST_QC.m does this)
%     2. dataDir/QC_result.mat written by computeCESTQualityMetrics, if
%        'UseSavedQC' is true (default)
%     3. none -> pre-QC behaviour only (as before), roiFit.qcApplied = false
%
%   FIGURES (per ROI, prefix = <SavePrefix><roiLabel>_)
%     if 'PlotMaps'      : <prefix>mask.png, <prefix>glucose_roi_overlay.png   (as before)
%     if 'PlotZspectrum' : <prefix>MeanROI_Zspectrum_with_LorentzFit.png       (pre-QC, as before)
%       and, when a QC mask is available:
%                          <prefix>MeanROI_Zspectrum_with_LorentzFit_postQC.png
%                          <prefix>MeanROI_Zspectrum_preVsPostQC.png
%
%   NAME-VALUE PAIRS
%     'RoiFile'          (default CESTQCConfig().RoiFile) - string OR cell array
%     'RotN'             (default 3)
%     'ClimBg'           (default [0, 150])
%     'ClimGl'           (default [0, 0.2])
%     'AlphaMask'        (default 0.3)
%     'AlphaGlucose'     (default 0.7)
%     'FigPosition'      (default [100 100 500 500])
%     'SavePrefix'       (default '')
%     'PlotMaps'         (default true)
%     'PlotZspectrum'    (default true)
%     'YLimGlucoseRight' (default [0, 0.2])
%     'BadPixelMask'     (default [])
%     'QCFlags'          (default [])
%     'UseSavedQC'       (default true)   load dataDir/QC_result.mat if no mask given
%     'SavePostQCFit'    (default true)   write dataDir/<SavePrefix>LorentzianFitPostQC.mat
%
%   MATLAB R2018a compatible.

def = CESTQCConfig();

p = inputParser;
addParameter(p, 'RoiFile', def.RoiFile);
addParameter(p, 'RotN', 3);
addParameter(p, 'ClimBg', [0, 150]);
addParameter(p, 'ClimGl', [0, 0.2]);
addParameter(p, 'AlphaMask', 0.3);
addParameter(p, 'AlphaGlucose', 0.7);
addParameter(p, 'FigPosition', [100, 100, 500, 500]);
addParameter(p, 'SavePrefix', '');
addParameter(p, 'PlotMaps', true);
addParameter(p, 'PlotZspectrum', true);
addParameter(p, 'YLimGlucoseRight', [0, 0.2]);
addParameter(p, 'BadPixelMask', []);
addParameter(p, 'QCFlags', []);
addParameter(p, 'UseSavedQC', true);
addParameter(p, 'SavePostQCFit', true);
parse(p, varargin{:});
opt = p.Results;

%% -------- load fit result --------
fitResultsPath = fullfile(dataDir, 'LorentzianFitResult.mat');
if ~isfile(fitResultsPath)
    error('plotCESTGlucoseROIResults:missingFit', ...
        'LorentzianFitResult.mat not found in %s. Run runCESTGlucoseROIFit first.', dataDir);
end
F = load(fitResultsPath, 'popt', 'P', 'gluc_integral_fullMap', 'M0_stack', 'Segment', 'Z_corrExt');
gluc_integral_fullMap = F.gluc_integral_fullMap;
M0_stack = F.M0_stack;
Z_corrExt = F.Z_corrExt;
W = F.P.SEQ.w;
[nx, ny] = size(gluc_integral_fullMap);

%% -------- QC mask source --------
badPixelMask = opt.BadPixelMask;
qcFlags = opt.QCFlags;
qcSource = 'none';
if ~isempty(badPixelMask)
    qcSource = 'argument';
elseif opt.UseSavedQC
    qcMatPath = fullfile(dataDir, [opt.SavePrefix 'QC_result.mat']);
    if isfile(qcMatPath)
        Q = load(qcMatPath, 'badPixelMask', 'qcFlags');
        badPixelMask = Q.badPixelMask;
        qcFlags = Q.qcFlags;
        qcSource = qcMatPath;
    end
end
qcApplied = ~isempty(badPixelMask);
if qcApplied
    badPixelMask = logical(badPixelMask);
    if ~isequal(size(badPixelMask), [nx, ny])
        error('plotCESTGlucoseROIResults:maskSizeMismatch', ...
            'badPixelMask is %s but gluc_integral_fullMap is %s in %s.', ...
            mat2str(size(badPixelMask)), mat2str([nx, ny]), dataDir);
    end
end
if ~isempty(qcFlags) && isfield(qcFlags, 'thresholdString')
    thresholdString = qcFlags.thresholdString;
elseif qcApplied
    thresholdString = 'thresholds: n/a (qcFlags not provided)';
else
    thresholdString = 'no QC applied';
end
fprintf('[plotCESTGlucoseROIResults] QC mask source: %s\n', qcSource);

%% -------- ROI list --------
if ischar(opt.RoiFile) || isstring(opt.RoiFile)
    roiFiles = {char(opt.RoiFile)};
else
    roiFiles = cellstr(opt.RoiFile);
end

nROI = numel(roiFiles);
roi_mean = nan(1, nROI);
roi_std = nan(1, nROI);
roi_labels = cell(1, nROI);
roiFit = struct('label', {}, 'file', {}, 'W', {}, 'qcApplied', {}, ...
    'nPix_pre', {}, 'meanZspec_pre', {}, 'popt_roi_pre', {}, 'gluc_integral_roiFit_pre', {}, ...
    'nPix_post', {}, 'meanZspec_post', {}, 'popt_roi_post', {}, 'gluc_integral_roiFit_post', {}, ...
    'roi_mean_pre', {}, 'roi_std_pre', {}, 'roi_mean_post', {}, 'roi_std_post', {});

M0_rot = rot90(M0_stack, opt.RotN);
gluc_rot = rot90(gluc_integral_fullMap, opt.RotN);

for k = 1:nROI
    roiFile = roiFiles{k};
    [~, roiLabel] = fileparts(roiFile);
    roiLabel = erase(roiLabel, 'ROI_');
    roi_labels{k} = roiLabel;
    prefix = [opt.SavePrefix roiLabel '_'];

    roiFit(k).label = roiLabel;
    roiFit(k).file = roiFile;
    roiFit(k).W = W;
    roiFit(k).qcApplied = qcApplied;
    roiFit(k).nPix_pre = NaN;  roiFit(k).meanZspec_pre = [];  roiFit(k).popt_roi_pre = [];  roiFit(k).gluc_integral_roiFit_pre = NaN;
    roiFit(k).nPix_post = NaN; roiFit(k).meanZspec_post = []; roiFit(k).popt_roi_post = []; roiFit(k).gluc_integral_roiFit_post = NaN;
    roiFit(k).roi_mean_pre = NaN; roiFit(k).roi_std_pre = NaN; roiFit(k).roi_mean_post = NaN; roiFit(k).roi_std_post = NaN;

    roiPath = fullfile(dataDir, roiFile);
    if ~isfile(roiPath)
        warning('plotCESTGlucoseROIResults:missingROI', 'ROI file not found, skipping: %s', roiPath);
        continue;
    end
    S = load(roiPath);
    if ~isfield(S, 'ROI')
        warning('plotCESTGlucoseROIResults:noROIvar', '%s has no variable ROI, skipping.', roiPath);
        continue;
    end
    ROI = S.ROI;
    if ~isequal(size(ROI), [nx, ny])
        error('plotCESTGlucoseROIResults:roiSizeMismatch', ...
            '%s is %s but gluc_integral_fullMap is %s in %s.', ...
            roiFile, mat2str(size(ROI)), mat2str([nx, ny]), dataDir);
    end
    ROI_logical = logical(ROI);
    if qcApplied
        ROI_clean = ROI_logical & ~badPixelMask;
    else
        ROI_clean = ROI_logical;
    end

    ROI_rot = rot90(ROI, opt.RotN);
    ROI_nan_rot = double(ROI_rot);
    ROI_nan_rot(ROI_nan_rot == 0) = NaN;
    croppedMap_rot = gluc_rot .* ROI_nan_rot;

    %% -------- ROI statistics (pixelwise map) --------
    roi_mean(k) = mean(croppedMap_rot(:), 'omitnan');
    roi_std(k) = std(croppedMap_rot(:), 'omitnan');
    roiFit(k).roi_mean_pre = roi_mean(k);
    roiFit(k).roi_std_pre = roi_std(k);
    if qcApplied
        cleanVals = gluc_integral_fullMap(ROI_clean);
        roiFit(k).roi_mean_post = mean(cleanVals(:), 'omitnan');
        roiFit(k).roi_std_post = std(cleanVals(:), 'omitnan');
        fprintf('[%s] pixel map  pre-QC: Mean %.4f, SD %.4f (n=%d) | post-QC: Mean %.4f, SD %.4f (n=%d)\n', ...
            roiLabel, roi_mean(k), roi_std(k), sum(ROI_logical(:)), ...
            roiFit(k).roi_mean_post, roiFit(k).roi_std_post, sum(ROI_clean(:)));
    else
        fprintf('[%s] Mean: %.4f, SD: %.4f\n', roiLabel, roi_mean(k), roi_std(k));
    end

    %% -------- Figures 1+2: mask overlay, glucose overlay (unchanged) --------
    if opt.PlotMaps
        fig1 = figure('Position', opt.FigPosition);
        imagesc(M0_rot, opt.ClimBg)
        colormap gray
        axis image
        hold on
        roi_red = cat(3, ones(size(ROI_rot)), zeros(size(ROI_rot)), zeros(size(ROI_rot)));
        h_roi = imshow(roi_red);
        set(h_roi, 'AlphaData', opt.AlphaMask * double(ROI_rot > 0));
        title(['mask - ' roiLabel], 'FontSize', 18);
        hold off
        set(gca, 'Position', [0.05, 0.05, 0.75, 0.9]);
        set(fig1, 'PaperPositionMode', 'auto');
        print(fig1, fullfile(dataDir, [prefix 'mask.png']), '-dpng', '-r300');

        fig2 = figure('Position', opt.FigPosition);
        ax1 = axes;
        axes(ax1);
        imagesc(M0_rot, opt.ClimBg)
        colormap(ax1, gray)
        axis image
        axis off
        ax2 = axes;
        axes(ax2);
        h_roi2 = imagesc(croppedMap_rot);
        colormap(ax2, jet)
        axis image
        set(h_roi2, 'AlphaData', opt.AlphaGlucose * double(~isnan(croppedMap_rot)));
        caxis(ax2, opt.ClimGl);
        cb = colorbar(ax2);
        cb.FontSize = 18;
        cb.Label.String = 'Glucose integral';
        ax2.XColor = 'none';
        ax2.YColor = 'none';
        ax2.Color = 'none';
        linkaxes([ax1, ax2])
        ax1.Position = ax2.Position;
        set(fig2, 'PaperPositionMode', 'auto');
        print(fig2, fullfile(dataDir, [prefix 'glucose_roi_overlay.png']), '-dpng', '-r300');
    end

    %% -------- mean ROI Z-spectrum fits: pre-QC and (if available) post-QC --------
    if opt.PlotZspectrum
        [meanZspec_pre, popt_roi_pre, P_pre, nPix_pre] = local_fitMeanZspec(Z_corrExt, ROI_logical, F.P);
        roiFit(k).nPix_pre = nPix_pre;
        roiFit(k).meanZspec_pre = meanZspec_pre;
        roiFit(k).popt_roi_pre = popt_roi_pre;
        roiFit(k).gluc_integral_roiFit_pre = local_glucIntegral(popt_roi_pre);

        % pre-QC figure, identical to the previous version
        local_plotPoolDecomposition(meanZspec_pre, popt_roi_pre, P_pre, W, opt, ...
            sprintf('Mean Z-spectrum - %s (pre-QC, n=%d)', roiLabel, nPix_pre), ...
            fullfile(dataDir, [prefix 'MeanROI_Zspectrum_with_LorentzFit.png']), '');

        if qcApplied
            nPix_post = sum(ROI_clean(:));
            if nPix_post == 0
                warning('plotCESTGlucoseROIResults:emptyPostQC', ...
                    '[%s] all %d ROI pixels were flagged by QC; post-QC fit skipped.', roiLabel, nPix_pre);
            else
                [meanZspec_post, popt_roi_post, P_post, nPix_post] = local_fitMeanZspec(Z_corrExt, ROI_clean, F.P);
                roiFit(k).nPix_post = nPix_post;
                roiFit(k).meanZspec_post = meanZspec_post;
                roiFit(k).popt_roi_post = popt_roi_post;
                roiFit(k).gluc_integral_roiFit_post = local_glucIntegral(popt_roi_post);

                fprintf('[%s] ROI-mean fit gluc integral  pre-QC: %.4f (n=%d) | post-QC: %.4f (n=%d)\n', ...
                    roiLabel, roiFit(k).gluc_integral_roiFit_pre, nPix_pre, ...
                    roiFit(k).gluc_integral_roiFit_post, nPix_post);

                % post-QC pool decomposition (same style as pre-QC figure)
                local_plotPoolDecomposition(meanZspec_post, popt_roi_post, P_post, W, opt, ...
                    sprintf('Mean Z-spectrum - %s (post-QC, n=%d)', roiLabel, nPix_post), ...
                    fullfile(dataDir, [prefix 'MeanROI_Zspectrum_with_LorentzFit_postQC.png']), thresholdString);

                % pre vs post comparison
                local_plotPreVsPost(meanZspec_pre, popt_roi_pre, P_pre, nPix_pre, ...
                    meanZspec_post, popt_roi_post, P_post, nPix_post, W, opt, roiLabel, ...
                    roiFit(k).gluc_integral_roiFit_pre, roiFit(k).gluc_integral_roiFit_post, ...
                    fullfile(dataDir, [prefix 'MeanROI_Zspectrum_preVsPostQC.png']), thresholdString);
            end
        end
    end
end

%% -------- save post-QC ROI fit --------
if opt.SavePostQCFit && qcApplied
    postQCPath = fullfile(dataDir, [opt.SavePrefix 'LorentzianFitPostQC.mat']);
    qcThresholds = [];
    if ~isempty(qcFlags) && isfield(qcFlags, 'thresholds')
        qcThresholds = qcFlags.thresholds;
    end
    generated = datestr(now);
    save(postQCPath, 'roiFit', 'badPixelMask', 'qcThresholds', 'thresholdString', 'qcSource', 'generated', '-v7');
    fprintf('Post-QC ROI fit saved: %s\n', postQCPath);
end

close all;

end


%% ===================== local functions =====================

function [meanZspec, popt_roi, P, nPix] = local_fitMeanZspec(Z_corrExt, mask, P)
% Mean Z-spectrum over mask (NaN-safe), then single Lorentzian multi-pool fit.
% P is taken fresh from LorentzianFitResult.mat by the caller because
% fitmodelfunc_NUM modifies it.
N_offsets = size(Z_corrExt, 4);
meanZspec = zeros(N_offsets, 1);
nPix = sum(mask(:));
for jj = 1:N_offsets
    slice_img = squeeze(Z_corrExt(:,:,1,jj));
    voxels = slice_img(mask);
    meanZspec(jj) = mean(voxels(~isnan(voxels)));
end
[P] = fitmodelfunc_NUM(meanZspec, P);
popt_roi = levmar_fit(meanZspec, P);
end


function g = local_glucIntegral(popt_roi)
% Same definition as gluc_integral_fullMap in runCESTGlucoseROIFit:
% amplitude (param 17) * width (param 18) * pi/2
g = popt_roi(17) * popt_roi(18) * pi/2;
end


function local_plotPoolDecomposition(meanZspec, popt_roi, P, W, opt, titleStr, pngPath, thresholdString)
[f, ~, f1, f2, f3, f4, f5, f6] = fitmodelfunc_ANA(popt_roi, P);
w_plot = linspace(min(W), max(W), 500)';

figZ = figure('Color', 'white', 'Position', [100 100 900 550]);
yyaxis left;
hold on;
plot(W, meanZspec, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot(w_plot, f(w_plot), '-k', 'LineWidth', 2);
plot(w_plot, f1(w_plot), '--c', 'LineWidth', 1.5);
plot(w_plot, f2(w_plot), '--b', 'LineWidth', 1.5);
plot(w_plot, f3(w_plot), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
plot(w_plot, f4(w_plot), '--g', 'LineWidth', 1.5);
plot(w_plot, f5(w_plot), '--r', 'LineWidth', 1.5);
ylabel('Z / Z_0', 'FontSize', 12, 'Color', 'k');
set(gca, 'YColor', 'k');
ylim([0 1]);

yyaxis right;
plot(w_plot, f6(w_plot), '-m', 'LineWidth', 2.5);
ylabel('Glucose pool amplitude', 'FontSize', 25, 'Color', 'm');
set(gca, 'YColor', 'm');
ylim(opt.YLimGlucoseRight);
hold off;

set(gca, 'FontSize', 25, 'XDir', 'reverse');
xlim([-4 4]);
xlabel('Offset (ppm)', 'FontSize', 25);
title(titleStr, 'FontSize', 16);
grid on;
box on;
if ~isempty(thresholdString)
    annotation(figZ, 'textbox', [0, 0.002, 1, 0.035], 'String', ['QC: ' thresholdString], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none', 'Interpreter', 'none', 'FontSize', 9);
end
set(figZ, 'PaperPositionMode', 'auto');
print(figZ, pngPath, '-dpng', '-r300');
end


function local_plotPreVsPost(meanZspec_pre, popt_roi_pre, P_pre, nPix_pre, ...
    meanZspec_post, popt_roi_post, P_post, nPix_post, W, opt, roiLabel, gPre, gPost, pngPath, thresholdString)
[f_pre, ~, ~, ~, ~, ~, ~, f6_pre] = fitmodelfunc_ANA(popt_roi_pre, P_pre);
[f_post, ~, ~, ~, ~, ~, ~, f6_post] = fitmodelfunc_ANA(popt_roi_post, P_post);
w_plot = linspace(min(W), max(W), 500)';
greyC = [0.55 0.55 0.55];

figC = figure('Color', 'white', 'Position', [100 100 900 550]);
yyaxis left;
hold on;
hDpre  = plot(W, meanZspec_pre, 'o', 'Color', greyC, 'MarkerSize', 5);
hFpre  = plot(w_plot, f_pre(w_plot), '--', 'Color', greyC, 'LineWidth', 1.5);
hDpost = plot(W, meanZspec_post, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
hFpost = plot(w_plot, f_post(w_plot), '-k', 'LineWidth', 2);
ylabel('Z / Z_0', 'FontSize', 12, 'Color', 'k');
set(gca, 'YColor', 'k');
ylim([0 1]);

yyaxis right;
hGpre  = plot(w_plot, f6_pre(w_plot), '--m', 'LineWidth', 1.5);
hGpost = plot(w_plot, f6_post(w_plot), '-m', 'LineWidth', 2.5);
ylabel('Glucose pool amplitude', 'FontSize', 20, 'Color', 'm');
set(gca, 'YColor', 'm');
ylim(opt.YLimGlucoseRight);
hold off;

set(gca, 'FontSize', 20, 'XDir', 'reverse');
xlim([-4 4]);
xlabel('Offset (ppm)', 'FontSize', 20);
title(sprintf('%s: mean ROI Z-spectrum, pre- vs post-QC\ngluc integral (ROI fit) pre %.4f (n=%d) | post %.4f (n=%d)', ...
    roiLabel, gPre, nPix_pre, gPost, nPix_post), 'FontSize', 14);
legend([hDpre, hFpre, hDpost, hFpost, hGpre, hGpost], ...
    {'data pre-QC', 'fit pre-QC', 'data post-QC', 'fit post-QC', 'glucose pool pre-QC', 'glucose pool post-QC'}, ...
    'Location', 'southwest', 'FontSize', 10);
grid on;
box on;
if ~isempty(thresholdString)
    annotation(figC, 'textbox', [0, 0.002, 1, 0.035], 'String', ['QC: ' thresholdString], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none', 'Interpreter', 'none', 'FontSize', 9);
end
set(figC, 'PaperPositionMode', 'auto');
print(figC, pngPath, '-dpng', '-r300');
end
