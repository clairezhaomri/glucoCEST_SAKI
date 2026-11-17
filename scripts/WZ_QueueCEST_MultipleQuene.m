%% MultipleQueue - 修正版
clear all;close all;
% clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');

%% 第一部分：處理 ctx ROI
fprintf('=== 開始處理 ctx ROI ===\n');
Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';

% 檢查路徑是否存在
if ~exist(Basepath, 'dir')
    error('Base path does not exist: %s', Basepath);
end

cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];

if isempty(SubjectPath)
    warning('No subjects found in %s', Basepath);
end

CESTROI.name=('ctx');

counter=0;
IntegralResultValue=cell(200,23);
IntegralResultHeader={'AnimalID','phase','diet','experimentType','AcquiredBaseOffsets','AcquiredBaseROIMeanValue','AcquiredInfuOffsets','AcquiredInfuROIMeanValue','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
     'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
     'rNOEBase-3.5ppm'};

for i=1:numel(SubjectPath)
    fprintf('Processing ctx - Subject %3.0f/%3.0f: %s\n', i, numel(SubjectPath), SubjectPath(i).name); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);
    
    if ~exist(SubjectMeasurementPath, 'dir')
        fprintf('Warning: Subject path does not exist: %s\n', SubjectMeasurementPath);
        continue;
    end
    
    cd(SubjectMeasurementPath);
    S = dir('*CEST_base*');
    
    for j=1:numel(S)         
        if S(j).isdir==1
            path=fullfile(S(j).folder,filesep,S(j).name);
            
            if ~exist(path, 'dir')
                fprintf('Warning: CEST path does not exist: %s\n', path);
                continue;
            end
            
            cd(path);
            counter=counter+1;
            
            % 提取實驗條件信息
            header=[S(j).folder];
            if contains(header,'Acute')
               IntegralResultValue{counter,2}='acute';
            elseif contains(header,'2wk')
               IntegralResultValue{counter,2}='2wk';
            elseif contains(header,'8wk')
               IntegralResultValue{counter,2}='8wk';
            end
            
            if contains(header,'Fat')
               IntegralResultValue{counter,3}='Fat';
            elseif contains(header,'Protein')
               IntegralResultValue{counter,3}='Protein';
            elseif contains(header,'Normal')
               IntegralResultValue{counter,3}='Normal';
            end

            if contains(header,'Feces')
               IntegralResultValue{counter,4}='Feces Infect';
            elseif contains(header,'Ctrl')
               IntegralResultValue{counter,4}='Ctrl Infect';
            end
            
            try
                name=CESTROI.name;
                [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,2,3.6,name);
                IntegralResultValue{counter,1}=AnimalID;
                IntegralResultValue{counter,9}=IntegralValue;
                
                [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,3,name);
                IntegralResultValue{counter,15}=Asym_SelectFreq;
                
                fprintf('  Processed: %s, Integral: %.4f, Asym: %.4f\n', AnimalID, IntegralValue, Asym_SelectFreq);
            catch ME
                fprintf('Error processing %s: %s\n', path, ME.message);
                IntegralResultValue{counter,1}='ERROR';
                IntegralResultValue{counter,9}=NaN;
                IntegralResultValue{counter,15}=NaN;
            end
        end
    end
end

% 保存ctx結果
cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTResult%s',date)));
if isfolder(ResultDataPath)==0
    mkdir(sprintf('glucoCESTResult%s',date));
end
cd(ResultDataPath);

save(sprintf('CESTAsymResultsGluco_%s_%s.mat',CESTROI.name,datestr(now,'yyyymmdd_HHMMSS')),'IntegralResultValue','IntegralResultHeader');
fprintf('=== ctx處理完成，共處理 %d 個樣本 ===\n\n', counter);

%% 第二部分：處理 medulla ROI
fprintf('=== 開始處理 medulla ROI ===\n');
clear IntegralResultValue; % 清除上一個ROI的結果

% 注意：這裡確認路徑是否正確
Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';

% 檢查路徑是否存在
if ~exist(Basepath, 'dir')
    error('Base path does not exist: %s', Basepath);
end

cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];

if isempty(SubjectPath)
    warning('No subjects found in %s', Basepath);
end

CESTROI.name=('medulla');

counter=0;
IntegralResultValue=cell(200,23);
IntegralResultHeader={'AnimalID','phase','diet','experimentType','AcquiredBaseOffsets','AcquiredBaseROIMeanValue','AcquiredInfuOffsets','AcquiredInfuROIMeanValue','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
     'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
     'rNOEBase-3.5ppm'};

for i=1:numel(SubjectPath)
    fprintf('Processing medulla - Subject %3.0f/%3.0f: %s\n', i, numel(SubjectPath), SubjectPath(i).name); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);
    
    if ~exist(SubjectMeasurementPath, 'dir')
        fprintf('Warning: Subject path does not exist: %s\n', SubjectMeasurementPath);
        continue;
    end
    
    cd(SubjectMeasurementPath);
    S = dir('*CEST_base*');
    
    for j=1:numel(S)         
        if S(j).isdir==1
            path=fullfile(S(j).folder,filesep,S(j).name);
            
            if ~exist(path, 'dir')
                fprintf('Warning: CEST path does not exist: %s\n', path);
                continue;
            end
            
            cd(path);
            counter=counter+1;
            
            % 提取實驗條件信息
            header=[S(j).folder];
            if contains(header,'Acute')
               IntegralResultValue{counter,2}='acute';
            elseif contains(header,'2wk')
               IntegralResultValue{counter,2}='2wk';
            elseif contains(header,'8wk')
               IntegralResultValue{counter,2}='8wk';
            end
            
            if contains(header,'Fat')
               IntegralResultValue{counter,3}='Fat';
            elseif contains(header,'Protein')
               IntegralResultValue{counter,3}='Protein';
            elseif contains(header,'Normal')
               IntegralResultValue{counter,3}='Normal';
            end

            if contains(header,'Feces')
               IntegralResultValue{counter,4}='Feces Infect';
            elseif contains(header,'Ctrl')
               IntegralResultValue{counter,4}='Ctrl Infect';
            end
            
            try
                name=CESTROI.name;
                [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,2,3.6,name);
                IntegralResultValue{counter,1}=AnimalID;
                IntegralResultValue{counter,12}=IntegralValue;
                
                [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,3,name);
                IntegralResultValue{counter,18}=Asym_SelectFreq;
                
                fprintf('  Processed: %s, Integral: %.4f, Asym: %.4f\n', AnimalID, IntegralValue, Asym_SelectFreq);
            catch ME
                fprintf('Error processing %s: %s\n', path, ME.message);
                IntegralResultValue{counter,1}='ERROR';
                IntegralResultValue{counter,12}=NaN;
                IntegralResultValue{counter,18}=NaN;
            end
        end
    end
end

% 保存medulla結果
cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTResult%s',date)));
if isfolder(ResultDataPath)==0
    mkdir(sprintf('glucoCESTResult%s',date));
end
cd(ResultDataPath);

save(sprintf('CESTAsymResultsGluco_%s_%s.mat',CESTROI.name,datestr(now,'yyyymmdd_HHMMSS')),'IntegralResultValue','IntegralResultHeader');
fprintf('=== medulla處理完成，共處理 %d 個樣本 ===\n', counter);

fprintf('\n=== 所有ROI批處理完成！ ===\n');