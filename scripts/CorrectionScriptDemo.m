addpath('/Users/clairezhao/Documents/MATLAB/CEST/');

E27=MRIDataLoader.LoadFrom2DSeq([pwd filesep 'pdata' filesep '1'])
E27.headerHistory{1}.method.MT_Offsets_NoM0

coregCEST=CESTdata;

coregCEST(:,:,3)=(coregCEST(:,:,1)+coregCEST(:,:,5))./2;

coregCEST(:,:,4)=(coregCEST(:,:,2)+coregCEST(:,:,6))./2;

coregCEST(:,:,5)=(coregCEST(:,:,3)+coregCEST(:,:,7))./2;
coregCEST(:,:,6)=(coregCEST(:,:,4)+coregCEST(:,:,8))./2;


coregCEST(:,:,7)=(coregCEST(:,:,5)+coregCEST(:,:,9))./2;
coregCEST(:,:,8)=(coregCEST(:,:,6)+coregCEST(:,:,10))./2;
coregCEST(:,:,9)=(coregCEST(:,:,7)+coregCEST(:,:,11))./2;
coregCEST(:,:,10)=(coregCEST(:,:,8)+coregCEST(:,:,12))./2;

coregCEST(:,:,11)=(coregCEST(:,:,9)+coregCEST(:,:,13))./2;
coregCEST(:,:,12)=(coregCEST(:,:,10)+coregCEST(:,:,14))./2;

coregCEST(:,:,13)=(coregCEST(:,:,11)+coregCEST(:,:,15))./2;
coregCEST(:,:,14)=(coregCEST(:,:,12)+coregCEST(:,:,16))./2;

coregCEST(:,:,15)=(coregCEST(:,:,13)+coregCEST(:,:,17))./2;

coregCEST(:,:,16)=(coregCEST(:,:,14)+coregCEST(:,:,18))./2;
% 
coregCEST(:,:,17)=(coregCEST(:,:,15)+coregCEST(:,:,19))./2;

 coregCEST(:,:,18)=(coregCEST(:,:,16)+coregCEST(:,:,20))./2;
coregCEST(:,:,19)=(coregCEST(:,:,17)+coregCEST(:,:,21))./2;

coregCEST(:,:,20)=(coregCEST(:,:,18)+coregCEST(:,:,22))./2;

coregCEST(:,:,21)=(coregCEST(:,:,19)+coregCEST(:,:,23))./2;
coregCEST(:,:,22)=(coregCEST(:,:,20)+coregCEST(:,:,24))./2;
% 
coregCEST(:,:,23)=(coregCEST(:,:,21)+coregCEST(:,:,25))./2;
coregCEST(:,:,24)=(coregCEST(:,:,22)+coregCEST(:,:,26))./2;
% 
coregCEST(:,:,25)=(coregCEST(:,:,23)+coregCEST(:,:,27))./2;
coregCEST(:,:,26)=(coregCEST(:,:,24)+coregCEST(:,:,28))./2;

% 
coregCEST(:,:,27)=(coregCEST(:,:,25)+coregCEST(:,:,29))./2;

coregCEST(:,:,28)=(coregCEST(:,:,26)+coregCEST(:,:,30))./2;
coregCEST(:,:,29)=(coregCEST(:,:,27)+coregCEST(:,:,31))./2;
 coregCEST(:,:,30)=(coregCEST(:,:,28)+coregCEST(:,:,32))./2;
  coregCEST(:,:,31)=(coregCEST(:,:,29)+coregCEST(:,:,33))./2;

coregCEST(:,:,32)=coregCEST(:,:,32);




ImageViewerStudio(coregCEST);
% 
% E8ca=E42c(:,:,:)./E42c(:,:,end);
% E19ca=E19(:,:,:)./E19(:,:,end);
% E42c=E8ca;
% E42c(:,:,11)=E19ca(:,:,1);
% E42c(:,:,12)=E19ca(:,:,2);

%% save correct data
save('CESTCorrectData','CESTCorrectData');

%%
save('coregCEST','coregCEST');
%%
path=pwd
%% Queue CESTMap
CESTAsymMap(path,0.5,2);

% %% Queue ROI_XXX.mat
%                  
% CESTROI.name=('ctx');
% ImageViewerStudio(CESTCorrectData);
% if isfile(sprintf('CESTROIValue_%s.mat',CESTROI.name))==1
%     load(sprintf('CESTROIValue_%s.mat',CESTROI.name))
%     preROI = ROIResult.roi(1).volume.vol;
%     ROI = rot90(preROI,1);
%     CESTROI.ROI=ROI;
%     save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
%     
% elseif isfile(sprintf('ROI_%s.mat',CESTROI.name))==0
%     preROI = ROIResult.roi(1).volume.vol;
%     ROI = rot90(preROI,1);
%     CESTROI.ROI=ROI;
%     save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
% end
%% Queue ROI_xxx.mat

CESTROI.name=('ctx');

% %for phantom
% CESTROI.name=('50mM');
% figure,imagesc(S0data);ROI=roipoly;
% CESTROI.ROI=ROI;
% save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
% %



if isfile(sprintf('CESTROIValue_%s.mat',CESTROI.name))==1
    load(sprintf('CESTROIValue_%s.mat',CESTROI.name))
    CESTROI.ROI=ROI;
    save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');

elseif isfile(sprintf('ROI_%s.mat',CESTROI.name))==0
    ROI=QueueROIDraw(path);
    CESTROI.ROI=ROI;
    % cd(path);
    save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
end
%% Queue Spectrum_xxx.mat
% load(sprintf('ROI_%s.mat',CESTROI.name));
if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
    if isfile(sprintf('ROI_%s.mat',CESTROI.name))==1
        load(sprintf('ROI_%s.mat',CESTROI.name));
        WZ_ROISpectrum(path,ROI,CESTROI.name);
    end
%     save(sprintf('spectrum_%s.mat',CESTROI.name),'-struct','CESTROI');
end


% clear all;


%% flirt commend
% flirt -in data.nii -ref sRef3D.nii -out rdata.nii -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -2D -dof 12  -interp trilinear
