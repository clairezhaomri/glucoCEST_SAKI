%% Queue - CEST fit + plot (base + infu), skips subjects already fit
clear all; close all;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
warning('off', 'curvefit:prepareFittingData:sizeMismatchInVariables');
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
if ~exist(Basepath, 'dir')
    error('Base path does not exist: %s', Basepath);
end

SubjectPath = [dir(fullfile(Basepath, '*MIBIO*')); dir(fullfile(Basepath, '*KHTCH*'))];
qcConfig = CESTQCConfig();
roiFiles = qcConfig.RoiFile;   % same ROI list as WZ_QueueCEST_QC.m

ResultHeader = {'AnimalID','session','phase','diet','experimentType','ROI','ROI_mean','ROI_std'};
ResultValue = cell(numel(SubjectPath) * 2 * numel(roiFiles) + 40, numel(ResultHeader));
counter = 0;

for i = 1:numel(SubjectPath)
    fprintf('Subject %3.0f/%3.0f: %s\n', i, numel(SubjectPath), SubjectPath(i).name);
    SubjectMeasurementPath = fullfile(SubjectPath(i).folder, SubjectPath(i).name);
    if ~exist(SubjectMeasurementPath, 'dir')
        continue;
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

        header = S(j).folder;
        phase = ''; diet = ''; experimentType = '';
        if contains(header, 'Acute'), phase = 'acute';
        elseif contains(header, '2wk'), phase = '2wk';
        elseif contains(header, '8wk'), phase = '8wk'; end
        if contains(header, 'Fat'), diet = 'Fat';
        elseif contains(header, 'Protein'), diet = 'Protein';
        elseif contains(header, 'Normal'), diet = 'Normal'; end
        if contains(header, 'Feces'), experimentType = 'Feces Infect';
        elseif contains(header, 'Ctrl'), experimentType = 'Ctrl Infect'; end

        fitPath = fullfile(path, 'LorentzianFitResult.mat');
        ok = true;
        try
            if isfile(fitPath)
                fprintf('  [%s/%s] already fit, skipping\n', SubjectPath(i).name, session);
            else
                runCESTGlucoseROIFit(path);
            end
        catch ME
            fprintf('  Fit error %s: %s\n', path, ME.message);
            ok = false;
        end

        if ~ok
            for r = 1:numel(roiFiles)
                counter = counter + 1;
                ResultValue(counter,:) = {SubjectPath(i).name, session, phase, diet, experimentType, ...
                    local_roi_label(roiFiles{r}), NaN, NaN};
            end
            continue;
        end

        try
            [roi_mean, roi_std, roi_labels] = plotCESTGlucoseROIResults(path, 'RoiFile', roiFiles);
        catch ME
            fprintf('  Plot error %s: %s\n', path, ME.message);
            roi_mean = nan(1, numel(roiFiles));
            roi_std = nan(1, numel(roiFiles));
            % BUGFIX: fileparts(f) with ONE output returns the folder (''),
            % not the name -> labels were always '' on the error path.
            roi_labels = cellfun(@(f) local_roi_label(f), roiFiles, 'UniformOutput', false);
        end

        for r = 1:numel(roi_labels)
            counter = counter + 1;
            ResultValue(counter,:) = {SubjectPath(i).name, session, phase, diet, experimentType, ...
                roi_labels{r}, roi_mean(r), roi_std(r)};
        end
    end
end

ResultValue = ResultValue(1:counter, :);

ReportBase = '/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/';
ResultDataPath = fullfile(ReportBase, sprintf('glucoCESTResult%s', date));
if ~isfolder(ResultDataPath), mkdir(ResultDataPath); end

save(fullfile(ResultDataPath, sprintf('CESTAnalyzeResults_%s.mat', datestr(now,'yyyymmdd_HHMMSS'))), ...
    'ResultValue', 'ResultHeader');
writetable(cell2table(ResultValue, 'VariableNames', ResultHeader), ...
    fullfile(ResultDataPath, sprintf('CESTAnalyzeResults_%s.csv', datestr(now,'yyyymmdd_HHMMSS'))));

fprintf('\n=== Analyze complete: %d ROI entries ===\n', counter);

function lbl = local_roi_label(f)
[~, name] = fileparts(f);
lbl = erase(name, 'ROI_');
end
