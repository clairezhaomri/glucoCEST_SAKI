clear all;close all;clc;

addpath('/Users/clairezhao/Documents/MATLAB/CEST')
Basepath='/Users/clairezhao/Downloads/series_19_fg_gre_FM_1.2x1.2x1.4_RR/';
cd(Basepath);
ProtocolPath=dir('MR.1.3*');
%% Rename DICOM file
 ContentTimeValue=cell(numel(ProtocolPath),2);
 counter=0;
for k=3:numel(ProtocolPath)
        fprintf('Iteration_counts:%3.3i/%3.3i\n', k, numel(ProtocolPath)); 

            FilePath=fullfile(Basepath,filesep,ProtocolPath(k).name);
            try%try to go into and load the data
            Data=MRIDataLoader.LoadFromDICOM(FilePath);
            NewName=Data.headerHistory{1}.contentTime;
            OldName=ProtocolPath(k).name;
            ContentTimeValue{k,1}=OldName;
            ContentTimeValue{k,2}=NewName;

%             movefile(OldName, sprintf('MR_%s',NewName));
            counter=counter+1;
            catch% if there is missing 2dseq...
                printf('No available data')
            end      
end
 


 
% 
% if S(16).name== '*Perfusion*'
%     disp('it works!')
% end
