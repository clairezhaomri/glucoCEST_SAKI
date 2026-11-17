clear all;close all;clc;

Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
% SubjectPath=dir('*KHTCH*');
%% Rename Bruker file
for i=1:numel(SubjectPath)
    
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path
    cd(SubjectMeasurementPath);
    S = dir('*CEST_infu*');%showing how many element are there in the folder     
     for j=1:numel(S)        
         
            try%try to go into and load the data
            ProtocolPath=fullfile(SubjectMeasurementPath,filesep,S(j).name);

            NewName='CEST_infu';
            OldPath=ProtocolPath;
            movefile(OldPath, sprintf('%s',NewName));
            catch% if there is missing 2dseq...
                printf('No available data')
            end   
     end

end
 

 
% %% test
% 
% if S(16).name== '*Perfusion*'
%     disp('it works!')
% end
