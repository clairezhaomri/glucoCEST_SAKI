clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/Reichenbach/MRI/PVL/';
cd(Basepath);
% SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
SubjectPath=dir('*TRANS*');



for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*Perfusion_FAIR_EPI_5*');%showing how many element are there in the folder     
     for j=1:numel(S)
            if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
                if isfile('Map_results.mat')==1
                    disp('Already Analyzed!')
                else
                    cd('/Users/clairezhao/brukermatlabscripts/FAIR_ASL_Perfusion')
                    try%try to go into and load the data                   
                        wz_ASL_T1Fitting(path)
                    catch% if there is missing 2dseq...                    
                        disp('No available data!')
                    end
                
                end
            end
     end
%      S = dir('*T1T2*');
%      for j=1:numel(S)
%             if S(j).isdir==1
%              path=fullfile(S(j).folder,filesep,S(j).name);
%              cd(path);
% % %                 if isfile('T1AcqAndLsqFit.mat')==1
% % %                      disp('Already Analyzed!')
% % %                 else
% % %                     cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit')
% % %                     try%try to go into and load the data                   
% % %                         T1MapFitting(path)
% % %                     catch% if there is missing 2dseq...                    
% % %                         disp('No available data!')
% % %                     end
% % %                 end
% %                 if isfile('T2AcqDataAndFit.mat')==1
% %                      disp('Already Analyzed!')
% %                 else
%                         cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit')
% %                         try%try to go into and load the data                   
%                         T2MapFitting(path)
% %                     catch% if there is missing 2dseq...                    
% %                         disp('No available data!')
% %                     end
%                 
% %                 end
%             end
%       end
     
    
end
%%
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];



for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
%      S = dir('*Perfusion*');%showing how many element are there in the folder     
%      for j=1:numel(S)
%             if S(j).isdir==1
%              path=fullfile(S(j).folder,filesep,S(j).name);
%              cd(path);
%                 if isfile('Map_results.mat')==1
%                     disp('Already Analyzed!')
%                 else
%                     cd('/Users/clairezhao/brukermatlabscripts/FAIR_ASL_Perfusion')
%                     try%try to go into and load the data                   
%                         wz_ASL_T1Fitting(path)
%                     catch% if there is missing 2dseq...                    
%                         disp('No available data!')
%                     end
%                 
%                 end
%             end
%      end
     S = dir('*T1T2*');
     for j=1:numel(S)
            if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%                 if isfile('T1AcqAndLsqFit.mat')==1
%                      disp('Already Analyzed!')
%                 else
%                     cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit')
%                     try%try to go into and load the data                   
%                         T1MapFitting(path)
%                     catch% if there is missing 2dseq...                    
%                         disp('No available data!')
%                     end
%                 end
                if isfile('T2AcqDataAndFit.mat')==1
                     disp('Already Analyzed!')
                else
                    cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit')
                    try%try to go into and load the data                   
                        T2MapFitting(path)
                    catch% if there is missing 2dseq...                    
                        disp('No available data!')
                    end
                
                end
            end
      end
     
    
end



%%
% if protocolPath.name=('*perfusion*')
%     disp('it works')
% end
% 
% 
% if isfile('Map_results.mat')==0
% 
%     disp('it works!')
% end
