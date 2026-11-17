%Wan-Ting Zhao. Universitätsklinikum Jena
%Version 2020. Oct.13

%% load data
clear all,clc;
% BasePath='/data/projects/bruker-hoerr-cest/data/ANEST-71670__HemeCEST_04_baseline__2019-10-02_20.28.59/';
% BasePath='/opt/PV6.0.1/data/mri/20201006_150311_MIBIO_00808_1_2/';


AnatomyPath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/Ctrl_Infections/20200319_144225_MIBIO_00546_1_2ProteinYoungCtrlnfektAcute/CEST_baseline_noTrig (E6)/';
AnatomyDataPath=fullfile(AnatomyPath,filesep,'pdata',filesep,'1');
AnatomyData=MRIDataLoader.LoadFrom2DSeq(AnatomyDataPath); 
AnatomySlice=squeeze(AnatomyData(:,:,end));
MTOffsets=AnatomyData.headerHistory{1}.method.MT_Offsets_NoM0;
cd(AnatomyPath);

if isfile ('BaseMap+-5ppm:spline')==1
    load('BaseMap+-5ppm:spline.mat');
       disp('cool')
elseif isfile('CESTBaseMap(spline).mat')==1
    load('CESTBaseMap(spline).mat')
end

% cd(BasePath);
% load('ROI_Kidney.mat');
% load('ROI_Phantom.mat')

ROIMeanAcqBase=zeros(size(MTOffsets,2),1);
ROIMeanBase=zeros(4000,1);
ROIstdAcqBase=zeros(size(MTOffsets,2),1);
ROIstdBase=zeros(4000,1);
[xSort,xIdx]=sort(MTOffsets);
MTsorted=MTOffsets(xIdx);
xaxisAcq=MTsorted;
xaxis=linspace(min(MTsorted),max(MTsorted),4000);%
%% ROI

%  CESTROIValue.name=('5ppmPhantom:12.5mM');
%  CESTROIValue.name=('5ppmMedulla');
 CESTROIValue.name=('medulla');
% 
cd (AnatomyPath);
if isfile(sprintf('ROI_%s.mat',CESTROIValue.name))==1
   load(sprintf('ROI_%s.mat',CESTROIValue.name));
   CESTROIValue.data=ROI;
else
    I=figure,imagesc(AnatomySlice);colormap gray;
    ROI1=roipoly;
    ROI2=roipoly;
    CESTROIValue.data=ROI1>0|ROI2>0;
    waitfor(I);
end



 for echoi = 1:length(Base_AcqData(:,1,1,1))
        %base trace
        singleechoBase=Base_AcqData(echoi,2,:,:);
        ROIMeanAcqBase(echoi)=mean2(singleechoBase(CESTROIValue.data));
        ROIstdAcqBase(echoi)=std2(singleechoBase(CESTROIValue.data));
 end



    for echoi = 1:length(Base_SplineCorrectedData(1,1,1,:))
        %base trace
        singleechoBase=Base_SplineCorrectedData(:,:,2,echoi);
        ROIMeanBase(echoi)=mean2(singleechoBase(CESTROIValue.data));
        ROIstdBase(echoi)=std2(singleechoBase(CESTROIValue.data));
%         VoxelNum=length(find(ROI_kidney>0));
%         ROISemBase=ROIstdBase./VoxelNum;
    end
    
%     %debug plot
%     figure,plot(xaxisAcq,ROIMeanAcqBase,'bo',xaxis,ROIMeanBase,'r-')
    
        %prepare data for aymmetric spectrum
        %make symmetric
       [min_value,min_idx]=min(ROIMeanBase);
         a=xaxis(1,min_idx);
         xaxisCorrected=xaxis-a;
        
        %generate symmetric datapoint
        synData_x_base=linspace(xaxisCorrected(1),xaxisCorrected(end),4000);
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
        

        
        %find idx
        start_Contrast = 200;
        startPositions_Contrast = find(synData_x_base<start_Contrast);
        startIdx = max(startPositions_Contrast)+1;
%                 
        stop_Contrast = 800;
        stopPositions_Contrast = find(synData_x_base<stop_Contrast);
        stopIdx = max(stopPositions_Contrast);
        
        Zero_Contrast = 0;
        ZeroPositions_Contrast = find(synData_x_base<Zero_Contrast);
        ZeroIdx = max(ZeroPositions_Contrast);
        
        xaxis4k=linspace(-2000,2000,4000);
        pp1=spline(synData_x_base',synData_y_base);
        pp2=spline(xaxis4k', synData_y1_base);

        y1=ppval(pp1,synData_x_base);%acquired data
        y2=ppval(pp2,synData_x_base);%synthesized symmetric data

%         asym_Base=abs((y2(startIdx:stopIdx)-(y1(startIdx:stopIdx))));%mirror trace-signal trace
%         asym_Base=abs(y2-y1);%mirror trace-signal trace
        asym_Base=abs(y2-y1);%mirror trace-signal trace
        asym_integral=trapz(asym_Base(startIdx:stopIdx));
%         %------------calcuate select frequency-----------------------------
%         GlucoCESTWorkHz = 1152;%264/523/832/1152
%         WorkPosition_Contrast_Positive = find(synData_x_base<GlucoCESTWorkHz);
%         WorkPositionIdx_Positive = max(WorkPosition_Contrast_Positive);
%         
%         WorkPosition_Contrast_Negative=find(synData_x_base<(-GlucoCESTWorkHz));
%         WorkPositionIdx_Negative = max(WorkPosition_Contrast_Negative)+1;
% 
%         %calculate Asymetric
%         glucoCESTSpectumValue_Positive=y1(WorkPositionIdx_Positive);
%         glucoCESTSpectumValue_Negative=y1(WorkPositionIdx_Negative);
%         
%         Asym_SelectFreq=glucoCESTSpectumValue_Negative-glucoCESTSpectumValue_Positive;
%         Asym_SelectFreq
%         %debug plot
%         figure,plot(synData_x_base,y1)
  
        %------------End of calcuating select frequency---------------------
        %save to structure
        CESTROIValue.ROIMeanBase=ROIMeanBase;
        CESTROIValue.ROIstdBase=ROIstdBase;
        CESTROIValue.ROIMeanAcqBase=ROIMeanAcqBase;
        CESTROIValue.ROIstdAcqBase=ROIstdAcqBase;
        CESTROIValue.xaxis1=xaxisCorrected;
        CESTROIValue.xaxis=xaxis;
        CESTROIValue.asym_Base=asym_Base;
        CESTROIValue.ROI=CESTROIValue.data;
        CESTROIValue.RawOffsets=MTOffsets';
        CESTROIValue.OffsetSorted=MTsorted';
        CESTROIValue.asym_integral=asym_integral;
        CESTROIValue.startContrast=start_Contrast;
        CESTROIValue.stopContrast=stop_Contrast;


cd (AnatomyPath);        
save(sprintf('CESTROIValue_%s.mat',CESTROIValue.name),'-struct','CESTROIValue');


%plot
  figure1=figure;
  axes1 = axes('Parent',figure1);
  
  plot1 = plot((xaxisAcq-a)./400,ROIMeanAcqBase,'b*','Linewidth',2); hold on;
  plot2 = plot((xaxisCorrected)./400, ROIMeanBase,'b-','Linewidth',2);
  YEqualZeroLine = zeros(length(xaxisCorrected),1);
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (ppm)','FontSize',15,'FontWeight','bold');
  
  
%   if stopIdx>4000 
  plot((synData_x_base(ZeroIdx:end)./400),asym_Base(ZeroIdx:end), 'm-','Linewidth',2);
  plot(xaxisCorrected, YEqualZeroLine,'k');
  
%   else
%   plot((synData_x_base(startIdx:stopIdx))./400,asym_Base, 'm-','Linewidth',2);
%   end
      
%   title(sprintf('ROI-based CEST & Asymmetric Spectrum(%s)',ROI.title),'FontSize',15);
  title(sprintf('CEST%s-%s',CESTData.header.patientID,CESTROIValue.name),'FontSize',15);
  legend({'AcqData','Spline','AsymSpectrum'},'FontSize',12, 'Location', 'northwest');
  legend('boxoff')
  
  box(axes1,'on');
  grid on;

%   
%   
%   if ROI.name='100ppm*'
%     axis([-100 100 0 1.2]);
%     set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
%     [-100 -50 -25 0 25 50 100],'YTick',...
%     [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2],'ZTick',[-1.2 0 1.2]);
%   elseif ROI.name=='12ppm*'
%      axis([-12 12 0 1]);
%     set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
%     [-12 -8 -4 0 4 8 12],'YTick',...
%     [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);
%   elseif ROI.name=='6ppm*'
    axis([-5 5 0 1]);
    set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [ 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2],'ZTick',[-1 0 1]);
%   elseif ROI.name=='3.5ppm*'
%     axis([-4 4 0 1.2]);
%     set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
%     [-4 -3 -2 -1 0 1 2 3 4],'YTick',...
%     [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);
%   end

saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string'))); 

%

ROIMeanAcqBaseCorrected=ROIMeanAcqBase;
save(sprintf('spectrum%s.mat',name),'ROIMeanAcqBaseCorrected','xaxisCorrected','synData_x_base','ROIMeanBase','asym_Base','asym_integral','xaxis');


   
   
   
% %% save figure
% % cd(BasePath);
%   saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
%   uisave({'xaxis_corrected1','xaxis_corrected2','ROIMeanBase','ROIMeanInfu','ROISemBase','ROISemInfu'},'ROI_Cortex_L_results');


     


