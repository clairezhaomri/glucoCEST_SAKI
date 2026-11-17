%% MultipleQueue
clear all;close all;
% clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];

% CESTROI.name=('ctx');
%%
% counter=0;
%      IntegralResultValue=cell(200,23);
%      IntegralResultHeader={'AnimalID','phase','diet','experimentType','AcquiredBaseOffsets','AcquiredBaseROIMeanValue','AcquiredInfuOffsets','AcquiredInfuROIMeanValue','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
%          'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
%          'rNOEBase-3.5ppm'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_infu*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%                   if isfolder('coregistration')
%                       cd coregistration
%                       delete coregCEST.mat
%                         if isfile('CESTCorrectData.mat')
%                         load('CESTCorrectData.mat')
%                         sRef3D=CESTCorrectData;
%                             for k=1:length(sRef3D(1,1,:))
%                                 sRef3D(:,:,k)=sRef3D(:,:,end);
%                             end
%                         MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
%                         save('sRef3D','sRef3D');
%                         else
%                         CESTdata=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']);
%                         sRef3D=CESTdata;
%                             for k=1:length(sRef3D(1,1,:))
%                                 sRef3D(:,:,k)=sRef3D(:,:,end);
%                             end
%                         MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
%                         save('sRef3D','sRef3D');
%                         end
%                         
%                         continue;
%                   else
%                       mkdir coregistration;
                      if isfile('CESTCorrectData.mat')
                        load('CESTCorrectData.mat')
                        cd coregistration
                        delete rdata.nii
                        MRIDataSaver.SaveToNIFTI('data.nii',CESTCorrectData)
                        sRef3D=CESTCorrectData;
                            for k=1:length(sRef3D(1,1,:))
                                sRef3D(:,:,k)=sRef3D(:,:,end);
                            end
                         MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
                         save('sRef3D','sRef3D');
                        
                      else
                        CESTdata=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']);
                        cd coregistration
                        delete rdata.nii
                        MRIDataSaver.SaveToNIFTI('data.nii',CESTdata)
                        sRef3D=CESTdata;
                            for k=1:length(sRef3D(1,1,:))
                                sRef3D(:,:,k)=sRef3D(:,:,end);
                            end
                        MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
                        save('sRef3D','sRef3D');
                        
                      end
           end
     end
 end







