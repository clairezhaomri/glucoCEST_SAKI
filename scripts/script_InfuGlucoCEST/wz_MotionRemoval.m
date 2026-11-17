%this script is prepared for removing motion slice out from CEST spectrum
%MotionCorrectedData=RAWDataSorted;
%MotionCorrectedData= 128x64x62;

%RAWOffset_Sorted= 1x62;
%6 11 14 19 24 25 34 42 47 50 54 60 to remove


RAWData=MRIDataLoader.LoadFrom2DSeq('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200318_084609_MIBIO_00559_1_1FatCtrlYoungAnimal/8CESTbase/pdata/1');
RAWOffset=RAWData.headerHistory{1}.method.MT_Offsets_NoM0;
RAWOffsetInfu=InfuData.headerHistory{1}.method.MT_Offsets_NoM0;
% RAWDataNormalized=RAWData(:,:,:)./RAWData(:,:,end);
% RAWDataNormalized(:,:,end)=[];



[xSort,xIndex]=sort(RAWOffset);
RAWOffset_Sorted=RAWOffset(xIndex(:));
RAWDataSorted=RAWData(:,:,xIndex);

[xSort,xIndex]=sort(RAWOffsetInfu);
RAWOffset_SortedInfu=RAWOffsetInfu(xIndex(:));
InfuDataSorted=InfuData(:,:,xIndex);
MotionCorrectedInfuData=InfuDataSorted;



ImageViewerStudio(RAWDataSorted);
MotionCorrectedData= RAWDataSorted;

% ToRemove=[6 11 14 19 24 25 34 42 47 50 54 60];%MIBIO-00549
ToRemove1=[3 5 7 10 11 27];%MIBIO-00486Base
ToRemove2=[1 3 15 20 22 23 27];%MIBIO-00486Infu


MotionCorrectedData(:,:,ToRemove(:))=[];
RAWOffset_SortedNew=RAWOffset_Sorted;
RAWOffset_SortedNew(ToRemove(:))=[];

MotionCorrectedInfuData(:,:,ToRemove2(:))=[];
RAWOffset_SortedInfuNew=RAWOffset_SortedInfu;
RAWOffset_SortedInfuNew(ToRemove2(:))=[];
ROIMeanAcqInfuMotionRemoval=ROIMeanAcqInfu;
ROIMeanAcqInfuMotionRemoval(ToRemove2(:))=[];


% ROIMeanAcqBaseCorrected=zeros(1,length(MotionCorrectedData(1,1,:)));
%  for echoi = 1:length(MotionCorrectedData(1,1,:))
%         %base trace
%         singleechoBase=MotionCorrectedData(:,:,echoi);
%         ROIMeanAcqBaseCorrected(echoi)=mean2(singleechoBase(ROI_kidneyCortex));
%         ROIMeanAcqBaseCorrected(echoi)=std2(singleechoBase(ROI_kidneyCortex));
%  end

%  
%  ROIMeanAcqBase=zeros(62,1);
%  for echoi = 1:length(Base_AcqData(:,1,1,1))
%         %base trace
%         singleechoBase=Base_AcqData(echoi,2,:,:);
%         ROIMeanAcqBase(echoi)=mean2(singleechoBase(ROI.data));
%         ROIstdAcqBase(echoi)=std2(singleechoBase(ROI.data));
%  end


 ROIMeanAcqBaseMotionRemoval=ROIMeanAcqBase;
 ROIMeanAcqBaseMotionRemoval(ToRemove(:))=[];
 


 xaxis=linspace(-1500,1500,100);
 BaseSpline=pchip(RAWOffset_SortedNew,ROIMeanAcqBaseMotionRemoval,xaxis)
 InfuSpline=pchip(RAWOffset_SortedInfuNew,ROIMeanAcqInfuMotionRemoval,xaxis);
 figure,plot(RAWOffset_SortedNew,ROIMeanAcqBaseMotionRemoval,'b*',xaxis,BaseSpline,'b-');hold on;
 plot(RAWOffset_SortedInfuNew, ROIMeanAcqInfuMotionRemoval,'r*',xaxis,InfuSpline,'r-');hold off;
 
 %%
 
       [min_value,min_idx]=min(BaseSpline);
         a=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-a;
        
        %generate symmetric datapoint
        synData_x_base=linspace(-1500,1500,100);
        synData_y_base=BaseSpline;

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
        stop_Contrast = 1500;
        stopPositions_Contrast1 = find(synData_x_base<stop_Contrast);
        stopIdx1 = max(stopPositions_Contrast1)+1;
% %-----with Infu----
%         %-----------------
        
        min_value=[];
        min_idx=[];
        FillNr=[];
        %-----------------
        [min_value,min_idx]=min(InfuSpline);
         a=xaxis(1,min_idx);
         xaxis_corrected2=xaxis-a; 
        %generate symmetric datapoint
        synData_x_Infu=linspace(-1500,1500,100);
        synData_y_infu=InfuSpline;

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
        stop_Contrast = 1500;
        stopPositions_Contrast2 = find(synData_x_Infu<stop_Contrast);
        stopIdx2 = max(stopPositions_Contrast2)+1;
%-----with Infu----        
%        correct display
        asym_Base(1:startIdx1)=0;
        asym_Base(stopIdx1:end)=0;
        asym_Infu(1:startIdx2)=0;
        asym_Infu(stopIdx2:end)=0;

       

 figure,plot(RAWOffset_SortedNew,ROIMeanAcqBaseMotionRemoval,'b*',xaxis,BaseSpline,'b-');hold on;
 plot(RAWOffset_SortedInfuNew, ROIMeanAcqInfuMotionRemoval,'r*',xaxis,InfuSpline,'r-');hold off;
        
        
%plot
  figure1=figure;
  axes1 = axes('Parent',figure1);
  
  plot1 = plot(RAWOffset_SortedNew,ROIMeanAcqBaseMotionRemoval,'*', xaxis,BaseSpline,'b-','Linewidth',2); hold on;
  plot2 = plot(RAWOffset_SortedInfuNew,ROIMeanAcqInfuMotionRemoval,'*',xaxis, InfuSpline,'r-','Linewidth',2);
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (Hz)','FontSize',15,'FontWeight','bold');
  axis([-1500 1500 0 1]);
  
  if stopIdx1>100 || stopIdx2>100
  plot(xaxis_corrected1(startIdx1:end),asym_Base(startIdx1:end), 'm-','Linewidth',2);
  plot(xaxis_corrected2(startIdx1:end),asym_Infu(startIdx1:end), 'k-','Linewidth',2);  
  else
  plot(xaxis_corrected1(startIdx1:stopIdx1),asym_Base(startIdx1:stopIdx1), 'm-','Linewidth',2);
  plot(xaxis_corrected2(startIdx1:stopIdx1),asym_Infu(startIdx1:stopIdx1), 'k-','Linewidth',2);
  end
      
  title('ROI-based CEST & Asymmetric Spectrum(cortex):MotionCorrect','FontSize',15);
  legend({'BaseAcqPt','Baseline','InfuAcqPt','Infusion','AsymBaseline','AsymInfusion'},'FontSize',15, 'Location', 'southwest');
  legend('boxoff')
  
  box(axes1,'on');


  set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-1500 -1000 -500 0 500 1000 1500],'YTick',...
    [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);


saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string')));  

  