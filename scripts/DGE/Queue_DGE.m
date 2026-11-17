%% Batch processing for DGE analysis (pixelwised map)
clear all; close all;

addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);

SubjectPath = dir('*MIBIO*');

for i = 1:numel(SubjectPath)
    fprintf('Processing animal %3.3i/%3.3i\n', i, numel(SubjectPath));
    
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
    
    try
        % Option 1: Single offset at 1.28 ppm
        GlucoCEST_DGE_SingleOffset(pathBase, pathInfu, 1.28);
        
        % Option 2: Band integration 0.5-2 ppm
%         GlucoCEST_DGE_BandIntegration(pathBase, pathInfu, [0.5, 2]);
        
        % Option 3: Band integration 1-3 ppm
%         GlucoCEST_DGE_BandIntegration(pathBase, pathInfu, [1, 3]);
        
    catch ME
        warning('Failed for %s: %s', SubjectPath(i).name, ME.message);
    end
    
end

close all;
