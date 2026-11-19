%% Queue - CEST QC only (RMSE/R2/cleaned ROI stats), no re-fit
% Run this ANYTIME after WZ_QueueCEST_Analyze.m has produced
% LorentzianFitResult.mat files. Independent, skippable, rerunnable.
%
% THRESHOLDS: none are defined in this script. computeCESTQualityMetrics
% takes its defaults from CESTQCConfig.m, the same file that
% plotCESTQCBeforeAfter uses interactively. Edit thresholds there only.
%
% PER SESSION this script does, in order, all with the SAME badPixelMask/qcFlags:
%   1. computeCESTQualityMetrics   -> QC_summary.txt, QC_result.mat
%   2. plotCESTQCBeforeAfter        -> <roi>_QC_before_after.png
%   3. plotCESTGlucoseROIResults    -> pre- vs post-QC mean-ROI Z-spectrum fits,
%                                      <roi>_MeanROI_Zspectrum_preVsPostQC.png,
%                                      <roi>_MeanROI_Zspectrum_with_LorentzFit_postQC.png,
%                                      LorentzianFitPostQC.mat
%
% REQUIRES on the MATLAB path (/Users/clairezhao/Documents/MATLAB/CEST/):
%   CESTQCConfig.m, computeCESTQualityMetrics.m, plotCESTQCBeforeAfter.m,
%   plotCESTGlucoseROIResults.m
% Verify with:  which CESTQCConfig computeCESTQualityMetrics plotCESTQCBeforeAfter plotCESTGlucoseROIResults
clear all; close all;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');

%% -------- run switches --------
forceRerun = true;   % true  -> recompute QC even if QC_summary.txt exists
                      %          (no need to delete QC_summary.txt by hand)
makePlots  = true;    % false -> skip plotCESTQCBeforeAfter (faster batch)
fitPostQC  = true;    % false -> skip post-QC mean-ROI Z-spectrum fit (plotCESTGlucoseROIResults)

%% -------- paths / config --------
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
if ~exist(Basepath, 'dir')
    error('Base path does not exist: %s', Basepath);
end

qcConfig = CESTQCConfig();          % single source of truth for thresholds
roiFiles = qcConfig.RoiFile;        % {'ROI_ctx.mat','ROI_medulla.mat'}

SubjectPath = [dir(fullfile(Basepath, '*MIBIO*')); dir(fullfile(Basepath, '*KHTCH*'))];

% One row per AnimalID x session x ROI (long on keys, wide on metrics).
%   group         = cohort prefix of AnimalID (MIBIO | KHTCH)
%   session       = base | infu   (from folder name, used for base-vs-infu pivots)
%   sessionFolder = raw session directory name (traceability, repeated sessions)
%   stamp         = this batch run, identical for all rows of one CSV
ResultHeader = {'AnimalID','group','session','sessionFolder','ROI', ...
                'meanRMSE','meanR2','badFraction', ...
                'ROI_n','ROI_nClean','ROI_mean_raw','ROI_std_raw','ROI_mean_clean','ROI_std_clean', ...
                'ROIfit_gluc_integral_pre','ROIfit_gluc_integral_post','stamp'};
ResultValue = cell(numel(SubjectPath) * 2 * numel(roiFiles) + 40, numel(ResultHeader));
counter = 0;
stamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('QC thresholds (CESTQCConfig):\n');
disp(qcConfig);

for i = 1:numel(SubjectPath)
    fprintf('Subject %3.0f/%3.0f: %s\n', i, numel(SubjectPath), SubjectPath(i).name);
    SubjectMeasurementPath = fullfile(SubjectPath(i).folder, SubjectPath(i).name);
    if ~exist(SubjectMeasurementPath, 'dir')
        continue;
    end
    group = 'other';
    if contains(SubjectPath(i).name, 'MIBIO')
        group = 'MIBIO';
    elseif contains(SubjectPath(i).name, 'KHTCH')
        group = 'KHTCH';
    end

    S = [dir(fullfile(SubjectMeasurementPath, '*CEST_base*')); ...
         dir(fullfile(SubjectMeasurementPath, '*CEST_infu*'))];

    for j = 1:numel(S)
        if S(j).isdir ~= 1
            continue;
        end
        path = fullfile(S(j).folder, S(j).name);

        session = 'base';
        if contains(S(j).name, 'infu', 'IgnoreCase', true)
            session = 'infu';
        end

        if ~isfile(fullfile(path, 'LorentzianFitResult.mat'))
            fprintf('  [%s/%s] no fit results yet, skipping\n', SubjectPath(i).name, session);
            continue;
        end

        qcPath = fullfile(path, 'QC_summary.txt');
        if isfile(qcPath) && ~forceRerun
            fprintf('  [%s/%s] QC already done, skipping (set forceRerun = true to redo)\n', ...
                SubjectPath(i).name, session);
            continue;
        end

        badPixelMask = [];
        qcFlags = [];
        try
            [meanRMSE, meanR2, badFraction, roi_mean_clean, roi_std_clean, roi_labels, badPixelMask, qcFlags] = ...
                computeCESTQualityMetrics(path, 'RoiFile', roiFiles);
        catch ME
            fprintf('  QC error %s: %s\n', path, ME.message);
            for si = 1:numel(ME.stack)
                fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
            meanRMSE = NaN; meanR2 = NaN; badFraction = NaN;
            roi_mean_clean = nan(1, numel(roiFiles));
            roi_std_clean  = nan(1, numel(roiFiles));
            roi_labels = cellfun(@(f) local_roi_label(f), roiFiles, 'UniformOutput', false);
        end

        if makePlots && ~isempty(badPixelMask)
            try
                % Reuse the SAME masks/flags as the summary -> figure and
                % QC_summary.txt can never disagree.
                plotCESTQCBeforeAfter(path, 'RoiFile', roiFiles, ...
                    'BadPixelMask', badPixelMask, 'QCFlags', qcFlags);
            catch ME
                fprintf('  QC before/after plot error %s: %s\n', path, ME.message);
                for si = 1:numel(ME.stack)
                    fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
                end
            end
        end

        % Post-QC mean-ROI Z-spectrum fit (pre vs post comparison + LorentzianFitPostQC.mat).
        % PlotMaps=false: mask/overlay PNGs were already made by WZ_QueueCEST_Analyze.m.
        roiFit = [];
        if fitPostQC && ~isempty(badPixelMask)
            try
                [~, ~, ~, roiFit] = plotCESTGlucoseROIResults(path, 'RoiFile', roiFiles, ...
                    'BadPixelMask', badPixelMask, 'QCFlags', qcFlags, 'PlotMaps', false);
            catch ME
                fprintf('  Post-QC ROI fit error %s: %s\n', path, ME.message);
                for si = 1:numel(ME.stack)
                    fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
                end
            end
        end

        for r = 1:numel(roi_labels)
            counter = counter + 1;
            ROI_n = NaN; ROI_nClean = NaN; ROI_mean_raw = NaN; ROI_std_raw = NaN;
            ROIfit_gluc_integral_pre = NaN; ROIfit_gluc_integral_post = NaN;
            if ~isempty(roiFit) && numel(roiFit) >= r
                ROIfit_gluc_integral_pre  = roiFit(r).gluc_integral_roiFit_pre;
                ROIfit_gluc_integral_post = roiFit(r).gluc_integral_roiFit_post;
            end
            if ~isempty(qcFlags) && isfield(qcFlags, 'roi') && numel(qcFlags.roi) >= r
                ROI_n        = qcFlags.roi(r).nROI;
                ROI_nClean   = qcFlags.roi(r).nClean;
                ROI_mean_raw = qcFlags.roi(r).meanRaw;
                ROI_std_raw  = qcFlags.roi(r).sdRaw;
            end
            ResultValue(counter,:) = {SubjectPath(i).name, group, session, S(j).name, roi_labels{r}, ...
                meanRMSE, meanR2, badFraction, ...
                ROI_n, ROI_nClean, ROI_mean_raw, ROI_std_raw, roi_mean_clean(r), roi_std_clean(r), ...
                ROIfit_gluc_integral_pre, ROIfit_gluc_integral_post, stamp};
        end
    end
end

ResultValue = ResultValue(1:counter, :);

%% -------- save batch report (results + thresholds used) --------
ReportBase = '/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/';
ResultDataPath = fullfile(ReportBase, sprintf('glucoCESTResult%s', date));
if ~isfolder(ResultDataPath), mkdir(ResultDataPath); end

save(fullfile(ResultDataPath, sprintf('CESTQCResults_%s.mat', stamp)), ...
    'ResultValue', 'ResultHeader', 'qcConfig', 'forceRerun');
writetable(cell2table(ResultValue, 'VariableNames', ResultHeader), ...
    fullfile(ResultDataPath, sprintf('CESTQCResults_%s.csv', stamp)));

% Thresholds used for this batch, as plain text next to the CSV
fid = fopen(fullfile(ResultDataPath, sprintf('CESTQCConfig_used_%s.txt', stamp)), 'w');
fprintf(fid, 'CESTQCConfig used by WZ_QueueCEST_QC on %s\n', datestr(now));
fprintf(fid, 'ConfigFile: %s\n', which('CESTQCConfig'));
fprintf(fid, 'forceRerun: %d | makePlots: %d | fitPostQC: %d\n\n', forceRerun, makePlots, fitPostQC);
fn = fieldnames(qcConfig);
for k = 1:numel(fn)
    v = qcConfig.(fn{k});
    if iscell(v)
        fprintf(fid, '%s: %s\n', fn{k}, strjoin(v, ', '));
    elseif islogical(v)
        fprintf(fid, '%s: %d\n', fn{k}, v);
    else
        fprintf(fid, '%s: %s\n', fn{k}, mat2str(v));
    end
end
fclose(fid);

fprintf('\n=== QC complete: %d ROI entries -> %s ===\n', counter, ResultDataPath);

function lbl = local_roi_label(f)
[~, name] = fileparts(f);
lbl = erase(name, 'ROI_');
end
