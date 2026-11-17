clear all,close all,clc;
%% load data
% BasePath='/data/projects/bruker-hoerr-cest/data/ANEST-71982__HemeCEST_06_Heme__2019-10-07_19.38.43/';
BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/ANEST_71983__HemeCEST_07_Heme__2019-10-08_05.07.44';
cd(BasePath);

T2hiresPath=dir('*T2RARE_highres_axialMS*');
T2hiresDataPath=fullfile(BasePath,T2hiresPath.name,filesep,'pdata',filesep,'1');
Hiresdata=MRIDataLoader.LoadFrom2DSeq(T2hiresDataPath); 
Hiresdata.Permute('read','line','echo');

load 'BaseMap.mat'
load 'InfuMap.mat'




% figure,imagesc(rot90(Hiresdata(:,:,8),3));
ROINum=1;
ROI=zeros(size(Base_SplineCorrectedData(:,1,1,1)),size(Base_SplineCorrectedData(1,:,1,1)),ROINum);
ROIMeanBase=zeros(1000,ROINum);
ROIMeanInfu=zeros(1000,ROINum);

ROIstdBase=zeros(1000,ROINum);
ROIstdInfu=zeros(1000,ROINum);

ROI(:,:,1) = tempROI;
for ROI_i=1:ROINum
    
    figure,imagesc(Base_SplineCorrectedData(:,:,1,1));
    ROI(:,:,ROI_i)=roipoly;
%         tempROI=roipoly;

%      ROI=tempROI(1:4:end,1:4:end);
%     ROI=tempROI(1:10:end,1:10:end);
%     ROI=rot90(ROI,1);%because imagesc has an flip image display,try to have same geometry hier.
    %for echoi = 1:length(Base_SplineCorrectedData(1,1,1,:))
        %base trace
        %singleechoBase=Base_SplineCorrectedData(:,:,2,echoi);
       % ROIMeanBase(echoi)=mean(singleechoBase(ROI));
        
        singleechoBase=squeeze(Base_SplineCorrectedData(:,:,2,:));
        ROI = (ROI).*1.0; % logical to double 
        ROI(ROI==0) = nan; % make zero numbers to nan
        ROI1 = repmat(ROI(:,:,ROI_i),[1 1 1000]);     
        ROIMeanBase(:,ROI_i) = squeeze(nanmean(nanmean(singleechoBase.*ROI)));
     
%         ROIstdBase(echoi)=std(singleechoBase(ROI));
%         %infusion trace
%         singleechoInfu=Infu_SplineCorrectedData(:,:,2,echoi);
%         ROIMeanInfu(echoi)=mean(singleechoInfu(ROI));
%         ROIstdInfu(echoi)=std(singleechoInfu(ROI));
    %end
    
end

figure,plot(xaxis, ROIMeanBase)

figure,imagesc(rot90(Hiresdata(:,:,8),3));
tempROI=roipoly;
ROI=tempROI(1:4:end,1:4:end);
ROI=rot90(ROI,1);%because imagesc has an flip image display,try to have same geometry hier.

