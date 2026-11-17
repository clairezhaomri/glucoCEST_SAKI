%% MultipleQueue
clear all;close all;

addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];

% CESTROI.name=('ctx');
%%
% for k=3:numel(SubjectPath)
for k=40:45
    fprintf('Iteration_counts:%3.3i/%3.3i\n', k, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(k).folder,SubjectPath(k).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_infu*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
             if isfolder('coregistration')
                      cd coregistration
%                         if isfile('coregCEST.mat')
%                             %cd(path);
%                             %run CEST asymmetric map with coregistered file
%                             %CESTAsymMap(path,0.5,2);
%                             % load exist ROI
%                             %load(sprintf('ROI_%s.mat',CESTROI.name));
%                             % get spectrum and save
%                             %WZ_ROISpectrum(path,ROI,CESTROI.name);
%                             continue;
%                         else
%                             cd coregistration
                            coregCEST=MRIDataLoader.LoadFromNIFTI('rdata.nii');
                            save('coregCEST','coregCEST');
                            cd(path);
                            %run CEST asymmetric map with coregistered file
                            CESTAsymMap(path,0.5,2);
                            % load exist ROI
                            CESTROI.name=('ctx');
                            load(sprintf('ROI_%s.mat',CESTROI.name));
                            % get spectrum and save
                            WZ_ROISpectrum(path,ROI,CESTROI.name);
                            CESTROI.name=('medulla');
                            load(sprintf('ROI_%s.mat',CESTROI.name));
                            % get spectrum and save
                            WZ_ROISpectrum(path,ROI,CESTROI.name);
                            
%                         end
                       
% 
%                       cd(path);
%                       cd ..
%                       S2=dir('*CEST_infu*');
%                       if numel(S2)==0  
%                           continue;
%                       else
%                       pathInfu=fullfile([S2.folder '/' S2.name]);
%                       cd(pathInfu);
%                         delete('WT-CESTasym*')
%                           if isfolder('coregistration')
%                                 continue;
%                           else
%                               mkdir coregistration;
%                               if isfile('CESTCorrectData.mat')
%                                 load('CESTCorrectData.mat')
%                                 S0=squeeze(CESTCorrectData(:,:,end));
%                                 Ref2D=S0;
%                                 Ref3D=repmat(Ref2D,[1,1,length(CESTCorrectData(1,1,:))]);
%                                 sRef3D=imgaussfilt(Ref3D,1);
%                                 cd coregistration
%                                 MRIDataSaver.SaveToNIFTI('data.nii',CESTCorrectData)
%                                 MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
%                                 save('Ref3D','Ref3D');
%                               else
%                                 CESTdata=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']);
%                                 S0=squeeze(CESTdata(:,:,end));
%                                 Ref2D=S0;
%                                 Ref3D=repmat(Ref2D,[1,1,length(CESTdata(1,1,:))]);
%                                 sRef3D=imgaussfilt(Ref3D,1);
%                                 cd coregistration
%                                 MRIDataSaver.SaveToNIFTI('data.nii',CESTdata)
%                                 MRIDataSaver.SaveToNIFTI('sRef3D.nii',sRef3D)
%                                 save('Ref3D','Ref3D');
%                               end
%                           end
             end
           end
                    
                  
              
      end
end




