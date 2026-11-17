%% MultipleQueue
clear all;close all;
% clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
% SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
 SubjectPath=dir('*MIBIO*');
%   SubjectPath=dir('*KHTCH*');


% CESTROI.name=('medulla');
%%
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_base*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
             try
%              SelectFreqAsymCESTMap(path,1.28);
             WZ_Replot_reverseXandTwoAxis(path);
             catch
             end
                     cd ..
                      S2=dir('*CEST_infu*');
                      pathInfu=fullfile([S2.folder '/' S2.name]);
                      cd(pathInfu);           
                       try
%                          SelectFreqAsymCESTMap(pathInfu,1.28);
                         WZ_Replot_reverseXandTwoAxis(pathInfu);
                       catch
                       end
           end
     end
end
%      close all;
