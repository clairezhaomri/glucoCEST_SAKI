function [Asym_SelectFreq] = WT_CESTSpectrumInSelectedHzMap(Datapath,glucoCESTppm,name,StartContrastPpm,EndContrastPpm)
%this function is designed to spline at selected freqeuncy band and store in  
CESTdata=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
        
cd(Datapath);
         
if isfile('CESTCorrectData.mat')
    load('CESTCorrectData.mat')
    CESTdata(:,:,:)=CESTCorrectData(:,:,:);
end

if isfolder('coregistration')
    cd('coregistration')
    load('coregCEST.mat')
    CESTdata(:,:,:)=coregCEST(:,:,:);
    cd ..
end

if isfile(sprintf('ROI_%s.mat',name))==1
    load(sprintf('ROI_%s.mat',name),'ROI'); 
    CESTdata=CESTdata.* ROI;
end

    CESTdata.Permute('echo','read','line');
    handles.CEST.CEST=CESTdata;
    handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;
    
    
    % check M0 location
if strcmp(CESTdata.headerHistory{1}.method.M0Option,'M0_Last')
    S0data=squeeze(CESTdata.vol(CESTdata.Size('echo'),:,:));
    CESTdata.vol(CESTdata.Size('echo'),:,:)=[];
    DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(CESTdata.Size('echo'));
    handles.CEST.DataSlope(CESTdata.Size('echo')) = [];
end

if strcmp(CESTdata.headerHistory{1}.method.M0Option,'M0_First')
    S0data=squeeze(CESTdata.vol(1,:,:));
    CESTdata.vol(1,:,:)=[];
    DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    handles.CEST.DataSlope(1) = [];
end
if strcmp(CESTdata.headerHistory{1}.method.M0Option,'NoM0')
    S0path='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/phantom_glucose_multiConcentr__CEST_multi_10__2019-05-07_14.23.36/10/';
    S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
    S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath); 
    S0data.Permute('echo','read','line');
    DataSlope_S0=S0data.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    S0data=squeeze(S0data(1,:,:));
end

    
 % make directory local
localCEST_MT_Offsets_NoM0 = handles.CEST.CEST.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localCEST_Vol = handles.CEST.CEST.vol;
localCEST_S0_Vol = S0data;
handles.CEST.DataSlopeS0 =DataSlope_S0;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


localSizeRead = handles.CEST.CEST.Size('read');
localSizeLine = handles.CEST.CEST.Size('line');   
    
AcqData = zeros(size(localCEST_MT_Offsets_NoM0,2),2,128,128);
SplineCorrectedData= zeros(128,128,2,4000);% splined x value,splined y value
SplineIntegralResults=zeros(128,128);%integral result    
AsymFreqResults=zeros(128,128);
AcqMTRasym=zeros(4000,128,128);
 for k = 1:localSizeRead
    fprintf('Iteration_Read %3.3i/%3.3i\n', k, localSizeRead); 
    for l = 1:localSizeLine
            if CESTdata(1,k,l)==0
            AcqData(:,:,k,l)=0;
            SplineCorrectedData(k,l,:,:)=0;
            SplineIntegralResults(k,l)=0;  
            AsymFreqResults(k,l)=0;  
            AcqMTRasym(:,k,l)=0;  
                continue;
            else
            %load data
            x = localCEST_MT_Offsets_NoM0;
               %this is for those with seperate measurement of S0 
%             y = ((localCEST_Vol(k,l,:) .* handles.CEST.DataSlope) ./ (localCEST_S0_Vol(k,l,1) .* handles.CEST.DataSlopeS0)).'; 
%             if localCEST_S0_Vol(k,l)==0
% %                 localCEST_S0_Vol(k,l)=(localCEST_S0_Vol(k-1,l));
%                 localCEST_S0_Vol(k,l)=(localCEST_S0_Vol(k,l));
%             else
            y = ((localCEST_Vol(:,k,l))  ./ (localCEST_S0_Vol(k,l)))'; %within same loop then DataSlope should be the same.
%             end   
          
            
            
            
            [xSort, xIndex] = sort(x);
            CEST_raw = [(x(xIndex(:))); y(xIndex(:))];
            CEST_corrected=zeros(size(CEST_raw));
            [min_val,min_idx]=min(CEST_raw(2,:));
              w0_FromMin=CEST_raw(1,min_idx);
%                 printf('Detected remaining minimum shift of CEST: %f ppm',w0_FromMin);
                %remove residual
             CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
             CEST_corrected(2,:) = CEST_raw(2,:);
      
                AcqData(:,1,k,l)= CEST_corrected(1,:);
                AcqData(:,2,k,l)= CEST_corrected(2,:);
                
                %work only for same step size within offset
%                 CEST_corrected_Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),100);% test effect of 100spline and 1000spline
               

                CEST_corrected_Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),4000);
                Spline_Y=spline(CEST_corrected(1,:)',CEST_corrected(2,:)',CEST_corrected_Spline_X);
      
      
                
                
        
        %find minimum as zero and correct
        [min_value,min_idx]=min(Spline_Y);
         B0Shift=CEST_corrected_Spline_X(1,min_idx);
         xaxis_corrected1=CEST_corrected_Spline_X-B0Shift;
            
        pp1=spline(xaxis_corrected1,Spline_Y);
        
        SplineCorrectedData(k,l,1,:)= xaxis_corrected1;
        SplineCorrectedData(k,l,2,:)= Spline_Y;
         
        xaxis=linspace(-2000,2000,4000);
        StartContrast1 = StartContrastPpm*400;
        startPositions_Contrast1 = find(xaxis<StartContrast1);
        startIdx = max(startPositions_Contrast1)+1;
                 
        EndContrast1 = EndContrastPpm*400;
        stopPositions_Contrast1 = find(xaxis<EndContrast1);
        stopIdx = max(stopPositions_Contrast1);

%         %debug plot        
%         figure, plot(xaxis,y1,'r',xaxis,y2,'b')
        Z=ppval(pp1,xaxis);

        %calculate asymSpectrum
        MTRasym=Z(end:-1:1)-Z;
%         asym_Base=(y2-y1);%mirror trace-signal trace     
%         IntegralValue=abs(trapz(asym_Base(startIdx:stopIdx)));
        xbin=zeros(1,4000);
        xbin(end/2+1:stopIdx)=1;
        IntegralValue=(trapz(xbin,MTRasym));    
        SplineIntegralResults(k,l)=IntegralValue; 

        
        
        
        %find the idx
        GlucoCESTWorkHz = glucoCESTppm*400;
        
        %calculate Asymetric
        glucoCESTSpectumValue_Positive=ppval(pp1,GlucoCESTWorkHz);
        glucoCESTSpectumValue_Negative=ppval(pp1,-GlucoCESTWorkHz);
        
        if glucoCESTppm==3.5%APT
            GlucoCESTWorkHz4ppm=4*400;
            GlucoCESTWorkHz3ppm=3*400;

            glucoCESTSpectumValue_4ppm=ppval(pp1,GlucoCESTWorkHz4ppm);
            glucoCESTSpectumValue_3ppm=ppval(pp1,GlucoCESTWorkHz3ppm);
            glucoCESTSpectumValue_APT=ppval(pp1,GlucoCESTWorkHz);
            Asym_SelectFreq=(glucoCESTSpectumValue_4ppm+glucoCESTSpectumValue_3ppm)./2-glucoCESTSpectumValue_APT;
            
        elseif glucoCESTppm==-3.5
            GlucoCESTWorkHzneg2ppm=-2*400;
            GlucoCESTWorkHz5ppm=5*400;
            
            glucoCESTSpectumValue_neg2ppm=ppval(pp1,GlucoCESTWorkHzneg2ppm);
            glucoCESTSpectumValue_5ppm=ppval(pp1,GlucoCESTWorkHz5ppm);
            glucoCESTSpectumValue_NOE=ppval(pp1,GlucoCESTWorkHz);
            Asym_SelectFreq=(glucoCESTSpectumValue_neg2ppm+glucoCESTSpectumValue_5ppm)./2-glucoCESTSpectumValue_NOE;
            
        else
               Asym_SelectFreq=glucoCESTSpectumValue_Negative-glucoCESTSpectumValue_Positive;
        end
        
    
   AsymFreqResults(k,l)=Asym_SelectFreq; 
   AcqMTRasym(:,k,l)=MTRasym;

   
   
   
            end
    end
% disp(fprintf('glucoCEST %3.2f ppm=%6.4f',glucoCESTppm, Asym_SelectFreq));
% save(sprintf('SelectFrqResult_%s.mat','glucoCESTppm'),'glucoCESTSpectumValue_Positive','glucoCESTSpectumValue_Negative' , 'Asym_SelectFreq','xaxis4k','ROIMeanBase4k');

% saveas(gcf,sprintf('CESTasym:1_2to2_4:%s_%s.png',datetime));
% savefig(gcf, sprintf('CESTasym:1_2to2_4%s.fig',char(datetime))); 
% cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/');
% ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
% 
% if isfolder(ResultDataPath)==0
% 
% mkdir(sprintf('glucoCESTppm_%3.2f',glucoCESTppm));
% else
% cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
% end
% 
% 
% ResultDataPathDate=fullfile(ResultDataPath,date);
% 
% 
% if isfolder(ResultDataPathDate)==0
%     mkdir(ResultDataPath, date);
%     cd([ResultDataPath,'/',date]);
% %     disp('0')
% else
%     cd([ResultDataPath,'/',date]);
% %     disp('1')
% end
% 
% 
% 
% save(sprintf('SelectFrqResult_%3.2fppm-%s-%s%s.mat',glucoCESTppm,Data.header.patientID, Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'GlucoCESTWorkHz', 'Asym_SelectFreq');

% save(sprintf('spectrum_1_2to2_4%d_%d_%s.mat',StartContrast,EndContrast,datetime),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis');

%additional attribute to daniel folder
% save(sprintf('spectrum_%s_%s_%s.mat',Data.header.patientID,Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis');


 end
 
    %Plot Mean Spectra 
   for echoi = 1:length(AcqData(:,1,1,1))
        %base trace
        singleecho=AcqData(echoi,2,:,:);
        ROIMeanAcq(echoi)=mean(singleecho(singleecho(:)~=0));
        ROIstdAcq(echoi)=std(singleecho(singleecho(:)~=0));
    end

    for echoi = 1:length(SplineCorrectedData(1,1,1,:))
        singleecho=SplineCorrectedData(:,:,2,echoi);
        ROIMeanSpline(echoi)=mean(singleecho(singleecho(:)~=0));
        ROIstdSpline(echoi)=std(singleecho(singleecho(:)~=0));
    end
    
     for echoi = 1:length(AcqMTRasym(:,1,1))
        singleecho=AcqMTRasym(echoi,:,:);
        MeanMTRasym(echoi)=mean(singleecho(singleecho(:)~=0));
        stdMTRasym(echoi)=std(singleecho(singleecho(:)~=0));
     end
     
%       pp2=spline(xaxis_corrected1,ROIMeanAcq);
%       meanZ=ppval(pp2,xaxis);
%       %calculate asymSpectrum
%       MeanMTRasym=meanZ(end:-1:1)-meanZ;
     
     residue=xaxis_corrected1(1)-CEST_corrected(1,1);
     
     figure1 = figure;
     % Create axes
     axes1 = axes('Parent',figure1);
     hold(axes1,'on');
     plot((CEST_corrected(1,:)+residue)./400, ROIMeanAcq,'*','Color',[0, 0.4470, 0.7410],'Linewidth',2);
%      plot((CEST_corrected(1,:)+residue)./400, ROIMeanAcq,'r*','Linewidth',2);

     hold on;
     plot(xaxis_corrected1./400,ROIMeanSpline,'-','Color',[0, 0.4470, 0.7410],'Linewidth',2);
%      plot(xaxis_corrected1./400,ROIMeanSpline,'r-','Linewidth',2);

     axis([-3.6 3.6 0 1.1])
     set(axes1,'FontSize',18,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [ 0  0.2  0.4  0.6  0.8  1.0],'ZTick',[-1 0 1]);

     ylabel('Water signal (S/S_0)','FontSize',15,'FontWeight','bold','Color',[0, 0.4470, 0.7410]);

     MeanMTRasym(stopIdx+1)=0;
     yyaxis right
     plot(xaxis(end/2+1:stopIdx+1)./400,MeanMTRasym(end/2+1:stopIdx+1),'-','Linewidth',2,'Color',[0.85 0.32 0.09])
%           plot(xaxis(end/2+1:stopIdx+1)./400,MeanMTRasym(end/2+1:stopIdx+1),'r-','Linewidth',2)

     ylim([0,0.65])
     ylabel('CEST signal','FontSize',18,'FontWeight','bold','Color',[0.85 0.32 0.09]);
     set(axes1,'YColor',[0.85 0.32 0.09],'FontSize',18,'FontWeight','bold');
     grid on;
     xlabel('Offset (ppm)','FontSize',18,'FontWeight','bold');
     title('Muiti-concentration glucose','FontSize',20);
     set(gca, 'XDir','reverse')
 
    
     %
     save('1mM.mat','residue','CEST_corrected','ROIMeanAcq','xaxis_corrected1','ROIMeanSpline','MeanMTRasym','stopIdx','xaxis','Asym_SelectFreq')
     
     
     
pause(10)

