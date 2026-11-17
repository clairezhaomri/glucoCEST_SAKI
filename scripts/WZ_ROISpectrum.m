%Wan-Ting Zhao. Universitätsklinikum Jena
%Version 2021. Jan.13
%This function provide the spectrum to be generated when there exist defined
%ROI as input(e.g.'ROI_ctx'/ 'ROI_ctx'). In general the pipeline will be as
%following:
%(a)queue and generate ROI,
%(b)queue and generate spectrum,
%(c)queue and generate multiple integral region/ selected Frequency, the
%entire dataset result will be store as handler and cells in workspace with
%timestamp


% load data
function [ROIMeanAcqBaseCorrected, ROIMeanBase, asym_Base,asym_integral, synData_x_base, xaxisAcq,xaxisCorrected]=WZ_ROISpectrum(path,ROI,name)
% clear all,clc;

DataPath=fullfile(path,filesep,'pdata',filesep,'1');
Data=MRIDataLoader.LoadFrom2DSeq(DataPath); 
MTOffsets=Data.headerHistory{1}.method.MT_Offsets_NoM0;

cd(path);
if isfile('BaseMap+-5ppm:spline.mat')==1
    load('BaseMap+-5ppm:spline.mat', 'Base_AcqData');
    load('BaseMap+-5ppm:spline.mat', 'Base_SplineCorrectedData');
else
    CESTAsymMap(path,0.5,2);
end


ROIMeanAcqBase=zeros(size(MTOffsets,2),1);
ROIMeanBase=zeros(4000,1);
ROIstdAcqBase=zeros(size(MTOffsets,2),1);
ROIstdBase=zeros(4000,1);
[xSort,xIdx]=sort(MTOffsets);
MTsorted=MTOffsets(xIdx);
xaxisAcq=MTsorted;
xaxis=linspace(min(MTsorted),max(MTsorted),4000);%

% ROI
 CESTROIValue.name=name;
% 
cd (path);
if isfile(sprintf('ROI_%s.mat',CESTROIValue.name))==1
   load(sprintf('ROI_%s.mat',CESTROIValue.name));
   CESTROIValue.data=ROI;
else
    disp('no proper ROI');
end



 for echoi = 1:length(Base_AcqData(:,1,1,1))
        %base trace
        singleechoBase=Base_AcqData(echoi,2,:,:);
        ROIMeanAcqBase(echoi)=mean2(singleechoBase(CESTROIValue.data));
        ROIstdAcqBase(echoi)=std2(singleechoBase(CESTROIValue.data));
 end

    for echoi = 1:length(Base_SplineCorrectedData(1,1,1,:))
        singleechoBase=Base_SplineCorrectedData(:,:,2,echoi);
        ROIMeanBase(echoi)=mean2(singleechoBase(CESTROIValue.data));
        ROIstdBase(echoi)=std2(singleechoBase(CESTROIValue.data));
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
        
         MTRasym=y1(end:-1:1)-y1;
         xbin=zeros(1,4000);
         xbin(startIdx:stopIdx)=1;
         asym_integral=(trapz(xbin,MTRasym)); 
        
        
%         asym_Base=abs(y2-y1);%mirror trace-signal trace
% %         asym_Base=(y2-y1);%mirror trace-signal trace
%         asym_integral=trapz(asym_Base(startIdx:stopIdx));
        
        
        
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

%update xaxisAcq
xaxisAcqCorrect=xaxisAcq-a;
        
cd(path);        
save(sprintf('CESTROIValue_%s.mat',CESTROIValue.name),'-struct','CESTROIValue');


%plot
  figure1=figure;
  axes1 = axes('Parent',figure1);
  
  plot1 = plot(xaxisAcqCorrect./400,ROIMeanAcqBase,'b*','Linewidth',2); hold on;
  plot2 = plot((xaxisCorrected)./400, ROIMeanBase,'b-','Linewidth',2);
  YEqualZeroLine = zeros(length(xaxisCorrected),1);
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (ppm)','FontSize',15,'FontWeight','bold');
  
  
  plot((synData_x_base(ZeroIdx:end)./400),asym_Base(ZeroIdx:end), 'm-','Linewidth',2);
  plot(xaxisCorrected, YEqualZeroLine,'k');

      
  title(sprintf('CEST%s-%s',Data.header.patientID,CESTROIValue.name),'FontSize',15);
  legend({'AcqData','Spline','AsymSpectrum'},'FontSize',12, 'Location', 'northwest');
  legend('boxoff')
  
  box(axes1,'on');
  grid on;


    axis([-5 5 0 1]);
    set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [ 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2],'ZTick',[-1 0 1]);


saveas(gcf, sprintf('new%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('new%s.fig',get(get(gca,'title'),'string'))); 

%

ROIMeanAcqBaseCorrected=ROIMeanAcqBase;
save(sprintf('spectrum_%s.mat',name),'ROIMeanAcqBaseCorrected','xaxisCorrected','synData_x_base','ROIMeanBase','YEqualZeroLine','ROIMeanAcqBase','ZeroIdx','asym_Base','asym_integral','xaxisAcqCorrect');
% delete BaseMap+-5ppm:spline.mat
close all;clear all;
end

