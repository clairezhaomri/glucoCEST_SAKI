%% load data
clear all,close all,clc;
% BasePath='/data/projects/bruker-hoerr-cest/data/ANEST-71670__HemeCEST_04_baseline__2019-10-02_20.28.59/';
BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/';
cd(BasePath);


% % 
% % T2hiresPath=dir('*T2RARE_highres_axialMS*');
% T2highresPath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/Untitled Folder/20191008_050744_ANEST_71983_1_1/5/';
% T2hiresDataPath=fullfile(T2highresPath,filesep,'pdata',filesep,'1');
% Hiresdata=MRIDataLoader.LoadFrom2DSeq(T2hiresDataPath); 
% Hiresdata.Permute('read','line','echo');

load 'BaseMap100.mat'
load 'InfuMap100.mat'
% load 'DiffMap.mat'



% figure,imagesc(rot90(Hiresdata(:,:,8),3));
% ROINum=1;
% ROI=zeros(length(Base_SplineCorrectedData(:,1,1,1)),length(Base_SplineCorrectedData(1,:,1,1)),ROINum);
ROIMeanBase=zeros(100,1);
ROIstdBase=zeros(100,1);

ROIMeanInfu=zeros(100,1);
ROIstdInfu=zeros(100,1);
xaxis=linspace(-1500,1500,100);%
%%

cd (BasePath);
if isfile('ROI.mat')
   load('ROI.mat');
else
    I=figure,imagesc(Base_SplineIntegralResults,[0 500]);colormap jet;
    ROI=roipoly;
    waitfor(I);
    save('ROI','ROI');
end

%      I=figure,imagesc(rot90(Hiresdata(:,:,8),3));
%        I=figure,imagesc(Hiresdata(:,:,8));
% 
%      tempROI=roipoly;
%      waitfor(I);
%      tempROI=tempROI(1:4:end,1:4:end);
%      tempROI=rot90(tempROI,1);%because imagesc has an flip image display,try to have same geometry hier.


    for echoi = 1:length(Base_SplineCorrectedData(1,1,1,:))
        %base trace
        singleechoBase=Base_SplineCorrectedData(:,:,2,echoi);
        ROIMeanBase(echoi)=mean2(singleechoBase(ROI));
        ROIstdBase(echoi)=std2(singleechoBase(ROI));
        VoxelNum=length(find(ROI>0));
        ROISemBase=ROIstdBase./VoxelNum;
        %infusion trace
        singleechoInfu=Infu_SplineCorrectedData(:,:,2,echoi);
        tempROIMeanInfu=singleechoInfu(ROI);
        ROIMeanInfu(echoi)=mean(mean(tempROIMeanInfu));
        ROIstdInfu(echoi)=std2(singleechoInfu(ROI));
        ROISemInfu=ROIstdInfu./VoxelNum; 
    end
    
    
    
%     uisave('tempROI','ROI_Medulla_L.mat');
        %make symmetric
       [min_value,min_idx]=min(ROIMeanBase);
         a=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-a;
        
        %generate symmetric datapoint
        synData_x_base=linspace(-1200,1200,120);
        synData_y_base=ROIMeanBase;

        synData_y1_base= synData_y_base;
        FillNr=(length(synData_y_base)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));

        %make sure 50pt in the left hand side
        if FillNr> 0
        synData_y1_base(1:FillNr-1)=synData_y_base(1:FillNr-1);
        synData_y1_base(FillNr:end)=synData_y_base(1:end-FillNr+1);% now symmetric
        synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip

       
        elseif FillNr< 0
        qFillNr=abs(FillNr);       
        synData_y1_base(1:qFillNr)=0;
        synData_y1_base(1:end-qFillNr)=synData_y_base(qFillNr+1:end);
        synData_y1_base(end-qFillNr+1:end)=synData_y_base(end-qFillNr+1:end);% now symmetric
        synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip
        
        elseif FillNr == 0
        synData_y1_base(1:end/2)=synData_y_base(1:end/2);
        synData_y1_base(end/2+1:end)=flipud(synData_y_base(1:end/2));
        end
        

        %calculate asymSpectrum
        asym_Base=abs(synData_y_base-synData_y1_base);
        
        %find display idx
        start_Contrast = 0;
        startPositions_Contrast1 = find(synData_x_base<start_Contrast);
        startIdx1 = max(startPositions_Contrast1)+1;
%                 
        stop_Contrast = 1800;
        stopPositions_Contrast1 = find(synData_x_base<stop_Contrast);
        stopIdx1 = max(stopPositions_Contrast1)+1;
% %-----with Infu----
%         %-----------------
        
        min_value=[];
        min_idx=[];
        FillNr=[];
        %-----------------
        [min_value,min_idx]=min(ROIMeanInfu);
         a=xaxis(1,min_idx);
         xaxis_corrected2=xaxis-a; 
        %generate symmetric datapoint
        synData_x_Infu=linspace(-1800,1800,100);
        synData_y_infu=ROIMeanInfu;

        synData_y1_Infu= synData_y_infu;
        FillNr=(length(synData_y_infu)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));

        %make sure 500pt in the left hand side
        if FillNr> 0
        synData_y1_Infu(1:FillNr-1)=synData_y_infu(1:FillNr-1);
        synData_y1_Infu(FillNr:end)=synData_y_infu(1:end-FillNr+1);% now symmetric
        synData_y1_Infu(end/2+1:end)=flipud(synData_y1_Infu(1:end/2));%now flip

       
        elseif FillNr< 0
        qFillNr=abs(FillNr);       
        synData_y1_Infu(1:qFillNr)=0;
        synData_y1_Infu(1:end-qFillNr)=synData_y1_Infu(qFillNr+1:end);
        synData_y1_Infu(end-qFillNr+1:end)=synData_y1_Infu(end-qFillNr+1:end);% now symmetric
        synData_y1_Infu(end/2+1:end)=flipud(synData_y1_Infu(1:end/2));%now flip
        
        elseif FillNr == 0
        synData_y1_Infu(1:end/2)=synData_y_infu(1:end/2);
        synData_y1_Infu(end/2+1:end)=flipud(synData_y_infu(1:end/2));
        end
        
        %calculate asymSpectrum
        asym_Infu=abs(synData_y_infu-synData_y1_Infu);
        %find display idx
        start_Contrast = 0;
        startPositions_Contrast2 = find(synData_x_Infu<start_Contrast);
        startIdx2 = max(startPositions_Contrast2)+1;
%                 
        stop_Contrast = 1800;
        stopPositions_Contrast2 = find(synData_x_Infu<stop_Contrast);
        stopIdx2 = max(stopPositions_Contrast2)+1;
%-----with Infu----        
%        correct display
        asym_Base(1:startIdx1)=0;
        asym_Base(stopIdx1:end)=0;
        asym_Infu(1:startIdx2)=0;
        asym_Infu(stopIdx2:end)=0;


%plot
  figure,errorbar(xaxis_corrected1, ROIMeanBase,ROISemBase,'r');hold on;   
  errorbar(xaxis_corrected2, ROIMeanInfu,ROISemInfu,'r');
  plot(xaxis_corrected1(startIdx1:stopIdx1),asym_Base(startIdx1:stopIdx1), 'b');
  plot(xaxis_corrected2(startIdx2:stopIdx2), asym_Infu(startIdx2:stopIdx2),'r');
%   plot(xaxis_corrected1(startIdx1:stopIdx1),asym_Base(startIdx1:stopIdx1), 'b--');
%   plot(xaxis_corrected2(startIdx2:stopIdx2), asym_Infu(startIdx2:stopIdx2),'r--');
  
  xlabel('Offset (Hz)');
  ylabel('S/S_0');
%   legend('Baseline','Infusion');
  axis([-1200 1200 0 1]);
%   title('sbj#71982_h_e_m_eCortex');
   title('vh170107:spline100');
   
%    saveas('tempROI',sprintf('ROI_%s.png',get(get,gca,'title'),'string'));

   

   
   
   
   
   
%% save figure
% cd(BasePath);
  saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
  uisave({'xaxis_corrected1','xaxis_corrected2','ROIMeanBase','ROIMeanInfu','ROISemBase','ROISemInfu'},'ROI_Cortex_L_results');


     



