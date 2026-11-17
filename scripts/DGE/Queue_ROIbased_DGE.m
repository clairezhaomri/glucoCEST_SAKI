%% Batch ROI-based DGE analysis
clear all; close all;

addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);

all_dirs = dir;
all_dirs = all_dirs([all_dirs.isdir]);


idx = contains({all_dirs.name}, 'Acute') & ...
      contains({all_dirs.name}, 'Normal') & ...
      contains({all_dirs.name}, 'Feces');
      
SubjectPath = all_dirs(idx);

% Define analysis parameters
SelectPpm = 1.28;  % Single offset
% SelectPpm = [0.5, 2];  % Or band integration

ROI_list = {'ROI_ctx', 'ROI_medulla'};  % Analyze both ROIs

% Collect results for group analysis
AllResults = struct();

for i = 1:numel(SubjectPath)
    fprintf('\n=== Processing animal %3.3i/%3.3i: %s ===\n', ...
        i, numel(SubjectPath), SubjectPath(i).name);
    
    SubjectMeasurementPath = fullfile(SubjectPath(i).folder, SubjectPath(i).name);
    cd(SubjectMeasurementPath);
    
    % Find baseline and infusion folders
    S_base = dir('*CEST_base*');
    S_infu = dir('*CEST_infu*');
    
    if isempty(S_base) || isempty(S_infu)
        warning('Missing baseline or infusion data for %s', SubjectPath(i).name);
        continue;
    end
    
    pathBase = fullfile(S_base(1).folder, S_base(1).name);
    pathInfu = fullfile(S_infu(1).folder, S_infu(1).name);
    
    % Process each ROI
    for r = 1:length(ROI_list)
        ROI_name = ROI_list{r};
        
        try
            fprintf('\n--- Analyzing %s ---\n', ROI_name);
            
            [AnimalID, DeltaS_Result] = GlucoCEST_DGE_ROI(pathBase, pathInfu, SelectPpm, ROI_name);
            
            % Store result
            if ~isfield(AllResults, ROI_name)
                AllResults.(ROI_name) = {};
            end
            AllResults.(ROI_name){end+1} = DeltaS_Result;
            
        catch ME
            warning('Failed for %s - %s: %s', SubjectPath(i).name, ROI_name, ME.message);
        end
    end
end

%% Group summary
fprintf('\n\n=== GROUP SUMMARY ===\n');
for r = 1:length(ROI_list)
    ROI_name = ROI_list{r};
    
    if isfield(AllResults, ROI_name) && ~isempty(AllResults.(ROI_name))
        fprintf('\n%s (n=%d):\n', ROI_name, length(AllResults.(ROI_name)));
        
        DeltaS_values = cellfun(@(x) x.DeltaS_percent, AllResults.(ROI_name));
        
        fprintf('  ΔS: %.2f ± %.2f%% (range: %.2f to %.2f%%)\n', ...
            mean(DeltaS_values), std(DeltaS_values), ...
            min(DeltaS_values), max(DeltaS_values));
    end
end

% Save all results
save([Basepath filesep 'AllResults_DGE_ROI.mat'], 'AllResults');

close all;
