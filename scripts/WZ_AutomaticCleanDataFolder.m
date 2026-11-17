%Wan-Ting Zhao; Universitaetsklinikum Jena; 2021.09.11
%This script is written for quality assurance. Basically it runs every folder to look for specific Term, if exist, the value would be loaded for some 
%quick calculated and document high standard deviated animals. The ROIs should be created under framework, which are represented as cell.   
%This script will generate a txt file for animals folders as a list. An entire list of T2 value would also be generated.

clear all;close all;clc;

Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/data_ready/';
cd(Basepath);
jointdir = [dir('*MIBIO*'); dir('*KHTCH*')];


%    counter=0;
%    DataBin=0;
   %TotalPerfusionValue =cell(numel(jointdir),4);%document a T2 ROI result list
%    fileID = fopen('AnimalList.txt','w');%create a txt file for motion data(Flirt required!)
%    fprintf(fileID,'%26s             %24s   %3s\n','path','AnimalName','Phase'); % print header of txt file
%% thresholding
for j=1:numel(jointdir)
    SubjectMeasurementPath=fullfile(jointdir(j).folder,jointdir(j).name);% measurement path    
    cd(SubjectMeasurementPath);
% 
     S1 = dir('CEST_base*');
     for k=1:numel(S1)
          if S1(k).isdir==1
           path=fullfile(S1(k).folder,filesep,S1(k).name);
           cd(path);
           
                  delete CESTROIValue_5ppmCtx.mat
                  delete BaseMap+-5ppm:spline.mat
          end
     end
%      end

     cd(SubjectMeasurementPath);             
     S3 = dir('*CEST_infu*');
%      counter=counter+1;
     for l=1:numel(S3)
          if S3(l).isdir==1
           pathInfu=fullfile(S3(l).folder,filesep,S3(l).name);
           cd(pathInfu);
                  delete CESTROIValue_5ppmCtx.mat
                  delete BaseMap+-5ppm:spline.mat
          end
     end
           
           
           
%             if isempty(S4)
%                 print(path);
%                 pause(5);
%                 fprintf(fileID,path);
                
%                 continue;
%                 
% %                 
%             else
%                 counter=counter+1;
%             origName = ('spectrum.mat');
%             newName  = ('spectrum_ctx.mat');
%             movefile(origName, newName);
          
end

% fclose(fileID);
cd(Basepath);
% TotalPerfusionValue(all(cellfun(@isempty, TotalPerfusionValue),2),:) = [];% remove empty cells
% save('Perfusion-medulla-valueList.mat','TotalPerfusionValue');
disp('Done!')



