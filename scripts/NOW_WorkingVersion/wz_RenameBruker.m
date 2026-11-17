clear all;close all;clc;

Basepath='/Volumes/T7Shield/MRIData/Joel_Data/';
cd(Basepath);
SubjectPath=dir('*ANEST*');
% SubjectPath=dir('*KHTCH*');
%% Rename Bruker file
for i=1:numel(SubjectPath)
    
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path
    cd(SubjectMeasurementPath);
    S = dir(SubjectMeasurementPath);%showing how many element are there in the folder

    for j=3:numel(S)%1st and 2nd element of folder are '.' and '..' for structure.
        if S(j).isdir==1 %if it is an directory
            ProtocolPath=fullfile(SubjectMeasurementPath,filesep,S(j).name);
            folderfullpath=fullfile(ProtocolPath,filesep,'pdata',filesep,'1');
            try%try to go into and load the data
            Data=MRIDataLoader.LoadFrom2DSeq(folderfullpath);
            NewName=Data.headerHistory{1}.protocolName;
            OldPath=ProtocolPath;
            NewNameReplacedCharacter=replace(NewName,"(","_");
            NewNameReplacedCharacter=replace(NewNameReplacedCharacter,")","");
            NewNameReplacedCharacter=replace(NewNameReplacedCharacter," ","");
            movefile(OldPath, sprintf('%s',NewNameReplacedCharacter));
            catch% if there is missing 2dseq...
                printf('No available data')
            end      
        end      
    end
end
 
disp('done!')
 
% %% test
% 
% if S(16).name== '*Perfusion*'
%     disp('it works!')
% end
