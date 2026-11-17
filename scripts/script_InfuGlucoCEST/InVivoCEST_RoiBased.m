%% load data
clear all;close all;clc;
% CEST part
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200127_115519_MIBIO_00484_1_2_protein/';
CESTpath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200127_115519_MIBIO_00484_1_2_protein/17_CEST31offsets/';
CESTdatapath=fullfile(CESTpath,filesep,'pdata',filesep,'1');
CESTdata=MRIDataLoader.LoadFrom2DSeq(CESTdatapath); 
CEST_Offsets=CESTdata.headerHistory{1}.method.MT_Offsets_NoM0;
CESTdata.Permute('read','line','echo');
handles.CEST.CEST=CESTdata;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


% %if no S0data
% S0=CESTdata.SubVolume('echo',1);
% DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(1);


%normalized_CESTdata=CESTdata./S0;
% % get separate S0 measurement 
% S0path='/data/projects/bruker-hoerr-cest/data/phantom_glucose_multiConcentr__CEST_multi_04__2019-02-06_17:53:37/wz_RARE_CEST_S0_E6';
% S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
% S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath);
% S0data.Permute('echo','read','line');
% handles.CEST.CEST_S0=S0data;


%% check M0 location
if strcmp(CESTdata.headerHistory{1}.method.M0Option,'M0_Last')
    S0data=CESTdata.vol(:,:,CESTdata.Size('echo'));
    CESTdata.vol(:,:,CESTdata.Size('echo'))=[];
    DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(CESTdata.Size('echo'));
    handles.CEST.DataSlope(CESTdata.Size('echo')) = [];
end

if strcmp(CESTdata.headerHistory{1}.method.M0Option,'M0_First')
    S0data=CESTdata.vol(:,:,1);
    CESTdata.vol(:,:,1)=[];
    DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    handles.CEST.DataSlope(1) = [];
end


%% make directory local
localCEST_MT_Offsets_NoM0 = handles.CEST.CEST.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localCEST_Vol = handles.CEST.CEST.vol;
localCEST_S0_Vol = S0data;
handles.CEST.DataSlopeS0 =DataSlope_S0;

localSizeRead = handles.CEST.CEST.Size('read');
localSizeLine = handles.CEST.CEST.Size('line');


%% manual draw Roi-based analysis
cd (Basepath);
if isfile('ROI.mat')
   load('ROI.mat');
elseif isfile('ROI_KidneyR.mat')
    load('ROI_KidneyR.mat');
else
    I=figure,imagesc(squeeze(CESTdata(:,:,1)),[0 30]);colormap gray;
    ROI=roipoly;
    waitfor(I);
    save('ROI','ROI');
end


[xSort, xIndex] = sort(CEST_Offsets);
locscaleddata=handles.CEST.CEST.vol ./ localCEST_S0_Vol;
locscaleddata=locscaleddata(:,:,xIndex);
    
RoiMean=zeros(handles.CEST.CEST.headerHistory{1}.method.Number_Offsets,1);
RoiStd=zeros(handles.CEST.CEST.headerHistory{1}.method.Number_Offsets,1);



for echoi = 1:handles.CEST.CEST.headerHistory{1}.method.Number_Offsets
   singleecho=locscaleddata(:,:,echoi);
   RoiMean(echoi)=mean(singleecho(ROI));
   RoiStd(echoi)=std(singleecho(ROI));
end

figure,plot(xSort,RoiMean,'b-o');
axis([-1500 1500 0 1]);
xlabel('frequency (Hz)');
ylabel('S/S_0');
title(sprintf('ROI-based Analysis:%s,%s',CESTdata.header.patientID,CESTdata.headerHistory{1}.visupars.VisuStudyId));
legend({'baseline'},'location','southwest');
% legend('baseline','infusion');
saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string')));

save('31OffsetROIResult.mat','RoiMean','xSort','ROI');
pause(1);



%% spline for ROIbased data
                

                cd(CESTpath);
                pp = spline(xSort(1,:)', R(2,:)');
                Spline_X=linspace(xSort(1,1),xSort(1,end),100);% test effect of 100spline and 1000spline
                Spline_Y=ppval(pp,Spline_X);
                [min_val,min_idx]=min(Spline_Y);
                w0_FromSplineMin=Spline_X(min_idx);
                Spline_Xcorrected=Spline_X - w0_FromSplineMin;
                %update AcqData x axis
                SplineData=zeros(size(Spline_X),2);

                SplineData(:,1)= xSort(:)-(w0_FromSplineMin./100);
                SplineData(:,2)= Spline_Y;
                
                start_Contrast = 400;
                startPositions_Contrast = find(Spline_X<start_Contrast);
                startIdx = max(startPositions_Contrast)+1;
%                 
                stop_Contrast = 1200;
                stopPositions_Contrast = find(Spline_X<stop_Contrast);
                stopIdx = max(stopPositions_Contrast)+1;
                
%                 if stopIdx > 100
%                     disp('exceed contrast');
%                     Base_SplineIntegralResults(k,l)=0;
%                     
%                 else
                    
                array2Trapz_x=Spline_X(startIdx:stopIdx);
                array2Trapz_y=1-(Spline_Y(startIdx:stopIdx));%flip upward for integration
                
                SplineIntegralResults=trapz(array2Trapz_x,array2Trapz_y);
                
%                 end
                
                
                
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
