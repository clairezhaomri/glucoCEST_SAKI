%Wan-Ting Zhao. Universitätsklinikum Jena
%Version 2020. Oct.13
function [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue]=WZ_QueueMedullaROIResult(path,StartContrastPpm,EndContrastPpm)

% load data
% clear all,clc; close all;

Data=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']); 
       ACQ_scan_name = Data.headerHistory{1}.acqp.ACQ_scan_name;
       ACQ_time = Data.headerHistory{1}.acqp.ACQ_time;
       AnimalID = Data.header.patientID;


cd(path)

if isfile('BaseMap+-5ppm:spline.mat')==1
   load('BaseMap+-5ppm:spline.mat');
    if isfile('Base_AcqData.mat')==0
% elseif isfile('CESTBaseMap(spline).mat')==1
%     load('CESTBaseMap(spline).mat')
% else
        cd('/Users/clairezhao/Documents/MATLAB/CEST')
        if size(Data.dims(1))==size(Data.dims(2))
            try
                CESTAsymMap(path,0.5,2)
            catch
                disp('dimension mismatch!')
            end
        end
    end
end


load('BaseMap+-5ppm:spline.mat');

if isfile('ROI_ctx.mat')==1
    load('ROI_ctx.mat');
else
    disp('no proper ROI!')
end

AnatomySlice=squeeze(Data(:,:,end));
MTOffsets=Data.headerHistory{1}.method.MT_Offsets_NoM0;


ROIMeanAcqBase=zeros(size(MTOffsets,2),1);
ROIMeanBase=zeros(4000,1);
ROIstdAcqBase=zeros(size(MTOffsets,2),1);
ROIstdBase=zeros(4000,1);
[xSort,xIdx]=sort(MTOffsets);
MTsorted=MTOffsets(xIdx);
xaxisAcq=MTsorted;
xaxis=linspace(min(MTsorted),max(MTsorted),4000);%


%ROI
CESTROI.name=('medulla');
CESTROI.data=ROI;




% extract ROIResult
 for echoi = 1:length(Base_AcqData(:,1,1,1))
        %base trace
        singleechoBase=Base_AcqData(echoi,2,:,:);
        ROIMeanAcqBase(echoi)=mean2(singleechoBase(CESTROI.data));
        ROIstdAcqBase(echoi)=std2(singleechoBase(CESTROI.data));
 end



    for echoi = 1:length(Base_SplineCorrectedData(1,1,1,:))
        %base trace
        singleechoBase=Base_SplineCorrectedData(:,:,2,echoi);
        ROIMeanBase(echoi)=mean2(singleechoBase(CESTROI.data));
        ROIstdBase(echoi)=std2(singleechoBase(CESTROI.data));
    end
    
       xaxis4k=linspace(-2000,2000,4000);%
       ROIMean4K=spline(ROIMeanAcqBase(1),ROIMeanAcqBase(end),xaxis4k);

    
%     %debug plot
%     figure,plot(xaxisAcq,ROIMeanAcqBase,'bo',xaxis,ROIMeanBase,'r-')
    
        %prepare data for aymmetric spectrum
        %make symmetric
       [min_value,min_idx]=min(ROIMean4K);
         B0Shift=xaxis4k(1,min_idx);
         xaxisCorrected=xaxis4k-B0Shift;
       RAW_x=xaxisCorrected;
       RAW_Y=ROIMeanAcqBaseCorrected;
        
        %generate symmetric datapoint
%         synData_x_base=linspace(-2000,2000,4000);
        synData_x_base=linspace(RAW_x(1),RAW_x(end),4000);

        synData_y_base=ROIMean4K;

        synData_y1_base= synData_y_base;
        FillNr=(length(synData_y_base)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));
        

        if isrow(synData_y_base)
            synData_y_base=synData_y_base';
        end
        
         if isrow(synData_y1_base)
            synData_y1_base=synData_y1_base';
        end

        %make sure half of the total pt in the left hand side
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
        


        %find integral interval of *mirror-ed* data
        StartContrast1 = StartContrastPpm*400;
        startPositions_Contrast1 = find(synData_x_base<StartContrast1);
        startIdx = max(startPositions_Contrast1)+1;
                 
        EndContrast1 = EndContrastPpm*400;
        stopPositions_Contrast1 = find(synData_x_base<EndContrast1);
        stopIdx = max(stopPositions_Contrast1);
        
        pp1=spline(synData_x_base',synData_y_base);
        pp2=spline(xaxis4k', synData_y1_base);

        y1=ppval(pp1,synData_x_base);
        y2=ppval(pp2,synData_x_base);
        
        
    
        
        %calculate asymSpectrum
        asym_Base=abs(y2-y1);%mirror trace-signal trace     
        IntegralValue=abs(trapz(asym_Base(startIdx:stopIdx)));
        
        CESTROI.ROIMeanBase=ROIMeanBase;
        CESTROI.ROIstdBase=ROIstdBase;
        CESTROI.ROIMeanAcqBase=ROIMeanAcqBase;
        CESTROI.ROIstdAcqBase=ROIstdAcqBase;
        CESTROI.xaxis1=xaxisCorrected;
        CESTROI.xaxis=xaxis;
        CESTROI.asym_Base=asym_Base;
        CESTROI.ROI=CESTROI.data;
        CESTROI.RawOffsets=MTOffsets';
        CESTROI.OffsetSorted=MTsorted';
        CESTROI.asym_integral=IntegralValue;
        CESTROI.startContrast=StartContrast1;
        CESTROI.stopContrast=EndContrast1;
        
        
        
        
        
        
        
        

% 
%     mkdir(Datapath, date);
%     cd([Datapath,'/',date]);
    
    
    disp(fprintf('Integral glucoCEST %3.1f-%3.1f ppm=%6.4f',StartContrastPpm,EndContrastPpm, IntegralValue))
%     save(sprintf('SelectFrqResult_%s.mat','glucoCESTppm'),'StartContrast','EndContrast' , 'IntegralValue','xaxis4k','ROIMeanBase4k');
    
%plot 
%   pause(10);
  h=figure; set(h,'Visible','off');
  axes1 = axes('Parent',h);
  
 
  plot((xaxisCorrected-B0Shift)./400,ROIMeanAcqBaseCorrected,'b*','Linewidth',2); hold on;
  plot((xaxisCorrected)./400, ROIMean4K,'b-','Linewidth',2)
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (ppm)','FontSize',15,'FontWeight','bold');
  
  
%   if stopIdx>4000
 plot((synData_x_base(end/2+1:end))./400,asym_Base(end/2+1:end), 'm-','Linewidth',2);   
%   else
%   plot((synData_x_base(startIdx:stopIdx))./400,asym_Base, 'm-','Linewidth',2);
%   end
      
%   title(sprintf('ROI-based CEST & Asymmetric Spectrum(%s)',ROI.title),'FontSize',15);
  title(sprintf('CESTMTRAsym:%s',AnimalID),'FontSize',15);
  legend({'AcqData','Spline','AsymSpectrum'},'FontSize',12, 'Location', 'southwest');
  legend('boxoff')
  
  box(axes1,'on');
  grid on;


    axis([-5 5 0 1]);
    set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);
% disp(fprintf('IntegralValue=%6.4f',IntegralValue))

saveas(gcf,sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.png',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name, ACQ_time));
savefig(gcf, sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.fig',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name,ACQ_time)); 
save(sprintf('spectrum:%3.1f-%3.1f_%s.mat',StartContrastPpm,EndContrastPpm,AnimalID,ACQ_scan_name, ACQ_time),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMean4K','asym_Base','IntegralValue','xaxis4k');




cd('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new');
ResultDataPath=fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new/',(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm)));

if isfolder(ResultDataPath)==0
mkdir(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm));
cd(ResultDataPath)
%     disp('0')
else
cd (fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new/',(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm))));
%     disp('1')
end


ResultDataPathDate=fullfile(ResultDataPath,date);


if isfolder(ResultDataPathDate)==0
    mkdir(ResultDataPath, date);
    cd([ResultDataPath,'/',date]);
%     disp('0')
else
    cd([ResultDataPath,'/',date]);
%     disp('1')
end


saveas(gcf,sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.png',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name, ACQ_time));
savefig(gcf, sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.fig',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name,ACQ_time)); 
save(sprintf('spectrum:%3.1f-%3.1f_%s%s%s.mat',StartContrastPpm,EndContrastPpm,AnimalID,ACQ_scan_name, ACQ_time),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMean4K','asym_Base','IntegralValue','xaxis4k');
close all;
end


     


