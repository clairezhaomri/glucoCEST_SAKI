%Wan-Ting Zhao, Universitätsklinikum Jena, Medizinische Mikrobiologie;
%version2020-05-02.This script is designed for Multiple Entry of CEST
%measurments, especially for abdominal imaging.User Must, (A). generate each
%measurement a CESTmap to continue. (B). offsets must incooperate
%at least 0(water frequency) and S0(6000 Hz) as internal control of each measurement. (C)On the
%occasion of corrupt S0,a seperate S0 and DataSlope are applicable.

%% load data
clear all,clc;
% BasePath='/data/projects/bruker-hoerr-cest/data/MIBIO-00539__NutritionCEST_69__2020-05-02_01.52.08/';


AnatomyPath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/20200107_120259_MIBIO_00482_1_1proteinAcute/CEST_baseline (E9)';
AnatomyDataPath=fullfile(AnatomyPath,filesep,'pdata',filesep,'1');
AnatomyData=MRIDataLoader.LoadFrom2DSeq(AnatomyDataPath); 
AnatomySlice=squeeze(AnatomyData(:,:,end));
MTOffsets=AnatomyData.headerHistory{1}.method.MT_Offsets_NoM0;
cd(AnatomyPath);


if isfile('CESTCorrectData.mat')
    load('CESTCorrectData.mat')
    AnatomyData=CESTCorrectData;
end

ROIMeanAcqBase=zeros(size(MTOffsets,2),1);
ROIMeanBase=zeros(4000,1);
ROIstdAcqBase=zeros(size(MTOffsets,2),1);
ROIstdBase=zeros(4000,1);

load ('BaseMap+-5ppm:spline.mat');
% load('CESTROIValue_5ppmCtx.mat');
% load('CESTROIValue_5ppmMedulla.mat');
[xSort,xIdx]=sort(MTOffsets);
MTsorted=MTOffsets(xIdx);
xaxis4k=linspace(min(MTsorted),max(MTsorted),4000);%
%% prepare correct data

ROIMeanAcqBaseCorrected=zeros(10,1);%depend on the size of real data
xaxisCorrected=zeros(size(ROIMeanAcqBaseCorrected));

%% when there is sufficiently good dataset

ROIMeanAcqBaseCorrected=ROIMeanAcqBase;
xaxisCorrected=OffsetSorted;

      %%  
        %work for multiple entry cest
        ROIMeanBase=spline(xaxisCorrected,ROIMeanAcqBaseCorrected,xaxis);
%         ROIMeanBase=spline(xaxisCorrected,ROIMeanAcqBaseCorrected,xaxis);

    
%     %debug plot
%     figure,plot(xaxisAcq,ROIMeanAcqBase,'bo',xaxis,ROIMeanBase,'r-')
    
        %prepare data for aymmetric spectrum
        %make symmetric
       [min_value,min_idx]=min(ROIMeanBase);
         a=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-a;
        
        %generate symmetric datapoint
        synData_x_base=linspace(xaxis_corrected1(1),xaxis_corrected1(end),4000);
        synData_y_base=ROIMeanBase;

        synData_y1_base= synData_y_base;
        FillNr=(length(synData_y_base)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));
        

        if isrow(synData_y_base)
            synData_y_base=synData_y_base';
        end
        
         if isrow(synData_y1_base)
            synData_y1_base=synData_y1_base';
        end

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
        


        
        %find display idx
        start_Contrast = 200;%1.2ppm
        startPositions_Contrast = find(synData_x_base<start_Contrast);
        startIdx = max(startPositions_Contrast)+1;
%                 
        stop_Contrast = 800;% 2.4 ppm
        stopPositions_Contrast = find(synData_x_base<stop_Contrast);
        stopIdx = max(stopPositions_Contrast);

        
        pp1=spline(synData_x_base',synData_y_base);
        pp2=spline(xaxis4k', synData_y1_base);

        y1=ppval(pp1,synData_x_base);
        y2=ppval(pp2,synData_x_base);

        
        
        
        %calculate asymSpectrum
        asym_Base=abs(y2-y1);%mirror trace-signal trace     
        asym_integral=(trapz(asym_Base(startIdx:stopIdx)));
        



%plot
  figure1=figure;
  axes1 = axes('Parent',figure1);
  
  plot1 = plot((xaxisCorrected-a)./400,ROIMeanAcqBaseCorrected,'b*','Linewidth',2); hold on;
  plot2 = plot((xaxis_corrected1)./400, ROIMeanBase,'b-','Linewidth',2);
  YAxisAtZero = zeros(length(xaxis_corrected1),1);
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (ppm)','FontSize',15,'FontWeight','bold');
  
  
%   if stopIdx>4000 
    plot((xaxis4k(end/2+1:end))./400,asym_Base(end/2+1:end), 'm-','Linewidth',2); 
    plot((xaxis_corrected1)./400,YAxisAtZero,'k-');

%   else
%   plot((xaxis(startIdx:stopIdx))./400,asym_Base, 'm-','Linewidth',2);
%   end
      
%   title(sprintf('ROI-based CEST & Asymmetric Spectrum(%s)',ROI.title),'FontSize',15);
  title(sprintf('correctedCEST%s%s',AnatomyData.header.patientID,name),'FontSize',15);
  legend({'AcqData','Spline','AsymSpectrum'},'FontSize',12, 'Location', 'southwest');
  legend('boxoff')
  
  box(axes1,'on');
  grid on;


     axis([-5 5 -0 1]);
    set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [-0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);


saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string'))); 
save(sprintf('spectrum%s.mat',name),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMeanBase','asym_Base','asym_integral','xaxis');
% cd ('/Users/clairezhao/brukermatlabscripts/CEST_DataAnalysis/Matlab_Analysis/InfusionGlucoCEST/');


   
   
   
% %% save figure
% % cd(BasePath);
%   saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
%   uisave({'xaxis_corrected1','xaxis_corrected2','ROIMeanBase','ROIMeanInfu','ROISemBase','ROISemInfu'},'ROI_Cortex_L_results');


     



