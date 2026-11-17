function [AnimalID,ACQ_scan_name,ACQ_time,Base_SplineIntegralResults] = CESTSelectFreqMap(path,StartContrastPpm,EndContrastPpm)
%Wan-Ting Zhao. Univerisitätsklinikum Jena
%Version 2020 Oct. 13
%This script spline and direct take minimum pixelwisely.

%        clear all;clc;
       CESTdata=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']);
       AnimalID = CESTdata.header.patientID;
       ACQ_scan_name = CESTdata.headerHistory{1}.acqp.ACQ_scan_name;
       ACQ_time = CESTdata.headerHistory{1}.acqp.ACQ_time;

if isfile('CESTCorrectData.mat')
    load('CESTCorrectData.mat')
    CESTdata(:,:,:)=CESTCorrectData(:,:,:);
end

if isfolder('coregistration')
    cd('coregistration')
    load('coregCEST.mat')
    CESTdata(:,:,:)=coregCEST(:,:,:);
end

    CESTdata.Permute('echo','read','line');
    handles.CEST.CEST=CESTdata;
    handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


% cd(BasePath);
% load('WASSRMap.mat');
% handles.WASSR.WASSRMap=WASSRMap; 





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
    S0path='/Volumes/ClaireData/BackUp20210108MacBook/Documents/data/MRI/pv6/CEST/RARECEST_PhantomValidateSeq/phantom_glucose_multiConcentr__CEST_multi_10__2019-05-07_14.23.36/10/';
    S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
    S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath); 
    S0data.Permute('echo','read','line');
    DataSlope_S0=S0data.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    S0data=squeeze(S0data(1,:,:));
end

% if 6k corrupt, also use external S0
% S0path='/opt/PV6.0.1/data/mri/20200502_015208_MIBIO_00539_1_2/16/';
%     S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
%     S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath); 
%     S0data.Permute('echo','read','line');
%     S0data=squeeze(S0data(1,:,:));
%     DataSlope_S0=S0data.headerHistory{1}.visupars.VisuCoreDataSlope(1);


% S0data = MRIDataBasic(S0data, 'read,line,slice,echo');







% make directory local
localCEST_MT_Offsets_NoM0 = handles.CEST.CEST.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localCEST_Vol = handles.CEST.CEST.vol;
localCEST_S0_Vol = S0data;
handles.CEST.DataSlopeS0 =DataSlope_S0;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


localSizeRead = handles.CEST.CEST.Size('read');
localSizeLine = handles.CEST.CEST.Size('line');

% %masking for faster looping
% cd (BasePath);
% if isfile('mask.mat')
%    load('mask.mat');
% else
%     mag_temp=squeeze(CESTdata(1,:,:));
%     I=figure,imagesc(mag_temp(:,:));
%     mask=roipoly;
%     waitfor(I);
%     save('mask','mask');
% end

 Base_AcqData = zeros(size(localCEST_MT_Offsets_NoM0,2),2,128,128);
 %TLfitData= zeros(2,1000,size(mask,1),size(mask,2));
 Base_SplineCorrectedData= zeros(128,128,2,4000);% splined x value,splined y value
 Base_SplineIntegralResults=zeros(128,128);%integral result
 %TLfitResults=zeros(size(mask,1),size(mask,2),1,3);
 
%pixelwise spline and store value  
tic
 for k = 1:localSizeRead
    fprintf('Iteration_Read %3.3i/%3.3i\n', k, localSizeRead); 
    for l = 1:localSizeLine
%             fprintf('Iteration_Line %3.3i/%3.3i\n', l, localSizeLine); 
%         if mask(k,l)== 0
%             handles.CEST.ContrastMatrixCorrected(k,l) = 0;
% %             handles.CEST.RMSEMap(k,l) = 0;
% %             handles.CEST.SSEMap(k,l) = 0;
% %            handles.CEST.residual(k,l) = 0;
%             Base_AcqData(:,:,k,l)=0;
%             %TLfitData(:,:,k,l)=0;
%             %TLfitResults(k,l,:,:)=0;
%             Base_SplineCorrectedData(k,l,:,:)=0;
%             Base_SplineIntegralResults(k,l)=0;
%             continue;
%         else
            x = localCEST_MT_Offsets_NoM0;
               %this is for those with seperate measurement of S0 
%             y = ((localCEST_Vol(k,l,:) .* handles.CEST.DataSlope) ./ (localCEST_S0_Vol(k,l,1) .* handles.CEST.DataSlopeS0)).'; 
            if localCEST_S0_Vol(k,l)==0
                localCEST_S0_Vol(k,l)=(localCEST_S0_Vol(k-1,l));
            else
            y = ((localCEST_Vol(:,k,l))  ./ (localCEST_S0_Vol(k,l)))'; %within same loop then DataSlope should be the same.
            end
            
                [xSort, xIndex] = sort(x);
                  CEST_raw = [(x(xIndex(:))); y(xIndex(:))];   % in Hz             
          
                %downsample point for testing
%                 x= downsample(x,5);
%                 y= downsample(y,5);

% -----------------------------direct take minimum ------------------------------------           
% % 
                  CEST_corrected=zeros(size(CEST_raw));
                  [min_val,min_idx]=min(CEST_raw(2,:));
                  w0_FromMin=CEST_raw(1,min_idx);
%                 printf('Detected remaining minimum shift of CEST: %f ppm',w0_FromMin);
                %remove residual
                  CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
                  CEST_corrected(2,:) = CEST_raw(2,:);
                  
% --------------------%for power calibration series ONLY --------------------------- 
%                    CEST_corrected(2,:) = CEST_corrected(2,:)./CEST_corrected(2,end);

%--------------------------------------Use WASSR----------------------------------
%                 CEST_corrected=zeros(size(CEST_raw));
%                 CEST_corrected(1,:) = CEST_raw(1,:) - handles.WASSR.WASSRMap(k,l);
%                 CEST_corrected(2,:) = CEST_raw(2,:);
                
% % %               %recheck residual after WASSR correction
% %                 [min_val,min_idx]=min(CEST_corrected(2,:));
% %                 w0_FromMin=CEST_corrected(1,min_idx);
% %                 printf('Detected remaining minimum shift of CEST: %f Hz',w0_FromMin);
% % %               %remove residual
% %                 CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
% %                 handles.CEST.residual(k,l) = w0_FromMin;
% %------------------------------------Single Lorentzian-----------------------------------------------------------------
%                 %prepare for initial values
% %                  [min_value,min_idx]=min(CEST_corrected(2,:));
% %                  w0=CEST_corrected(1,min_idx);
%                  Ampl0=mean(CEST_corrected(2,1:6));
%                  
%                 %lorentzian
%                  beta0=[Ampl0,0,0.03,1];%find minimum for initial values
%                  beta=nlinfit(CEST_corrected(1,1:min_idx+2),CEST_corrected(2,1:min_idx+2)',@lorentzian,beta0);
%                  Lfit= lorentzian(beta,xx);
%                  
%                  %goodness of fitting
% %                 square_error= (CEST_corrected(2,:)'-lorentzian(beta,CEST_corrected(1,:))).^ 2;
% %                 RMSE= sqrt(mean(square_error));
% %                 SSE= sqrt(sum(square_error));%sum of standard error
% %                 handles.CEST.RMSEMap(k,l) = RMSE;
% %                 handles.CEST.SSEMap(k,l) = SSE;
%                 
% %                 %find minimum
% %                  yMinValue = min(Lfit);
% %                  xPos = find(Lfit == yMinValue);
% %                  xValue = xx(xPos);
% %                  
%                 pp = spline(CEST_corrected(1,:)', CEST_corrected(2,:)');
% %                 zSpectra_corrected = [xx; lorentzian(beta,xx)'];
%                 qq = ppval(pp,xx)-lorentzian(beta,xx)';
%                 asymm_corrected =[xx;qq];
% %                 asymm_corrected(1,:) = xx;
% %                 asymm_corrected(2,:) = ppval(pp,xx)-lorentzian(beta,xx)';
% 
% 
% %                 for m = 1:1000
% %                     pos = lastValue * m/1000;
% %                     asymm_corrected(1,m) = pos;
% %                     asymm_corrected(2,m) = ppval(pp,-xx:-xPos)-lorentzian(beta,xx(xPos+1:end));
% %                 end
%                 
%                 %calculate goodness of fitting
% %                 square_error=((ppval(pp,xx(xPos-50:xPos+50)))'-lorentzian(beta,xx(1,xPos-50:xPos+50))).^ 2;
% %                 RMSE= sqrt(mean(square_error));
% %                 SSE= sqrt(sum(square_error));%sum of standard error
% %                 handles.CEST.RMSEMap(k,l) = RMSE;
% %                 handles.CEST.SSEMap(k,l) = SSE;
% %                 
%                 
%                 %store value
                Base_AcqData(:,1,k,l)= CEST_corrected(1,:);
                Base_AcqData(:,2,k,l)= CEST_corrected(2,:);
%                 
                
%------------------------------------Multiple Lorentzian-----------------------------------------------------------------
                                 
                %for power calibration series only
%                 [params,xFit,yFit,lorentzInt] = TripleLorentzFit(CEST_corrected(1,1:end-1),CEST_corrected(2,1:end-1));

%               %TripleLorentzFit  
%                 [params,xFit,yFit,lorentzInt] =CEST.TripleLorentzFit(CEST_corrected(1,:),CEST_corrected(2,:));
%                 TLfitResults(k,l,:,:)= [lorentzInt(1);lorentzInt(2);lorentzInt(3)];
%                 TLfitData(:,:,k,l)= [xFit;yFit];
%                 handles.CEST.ContrastMatrixCorrected(k,l,1,2)= TLfitResults(k,l,1,2)+TLfitResults(k,l,1,3);
                
%-----------------------------------------SPLINE------------------------------------------------------------------------- 
                CEST_corrected_X=(CEST_corrected(1,:)');
                CEST_corrected_Y=(CEST_corrected(2,:)');
                
                %work only for same step size within offset
%                 CEST_corrected_Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),100);% test effect of 100spline and 1000spline
               

                CEST_corrected_Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),4000);
%                 Spline_Y=pchip(CEST_corrected_X,CEST_corrected_Y,CEST_corrected_Spline_X);
                
                Spline_Y=spline(CEST_corrected_X,CEST_corrected_Y,CEST_corrected_Spline_X);

                
%                 %debug plot
%                 figure,plot(CEST_corrected_X,CEST_corrected_Y,'b*');
%                 hold on;
%                 plot(CEST_corrected_Spline_X, Spline_Y,'r-');
%                 hold off;

                [min_val,min_idx]=min(Spline_Y);
                w0_FromSplineMin=CEST_corrected_Spline_X(min_idx);
                Spline_Xcorrected=CEST_corrected_Spline_X - w0_FromSplineMin;
                Base_SplineCorrectedData(k,l,1,:)= Spline_Xcorrected;
                Base_SplineCorrectedData(k,l,2,:)= Spline_Y;
                
                %use spline data for symmetricity
                synData_x_base=Spline_Xcorrected;
                synData_y_base=Spline_Y;
                
               synData_y1_base= synData_y_base;
                FillNr=(length(synData_y_base)./2)-min_idx;
        

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
        


%                 %update AcqData x axis
                Base_AcqData(:,1,k,l)= CEST_corrected(1,:)'-(w0_FromSplineMin);
                
                    
                    start_Contrast = StartContrastPpm .*400;%1.2ppm
                    startPositions_Contrast = find(synData_x_base<start_Contrast);
                    startIdx = max(startPositions_Contrast)+1;
%                 
                    stop_Contrast = EndContrastPpm .*400;% 2.4 ppm
                    stopPositions_Contrast = find(synData_x_base<stop_Contrast);
                    stopIdx = max(stopPositions_Contrast);
                    xaxis4k=linspace(-2000,2000,4000);

                        
                    pp1=spline(synData_x_base,synData_y_base);
                    pp2=spline(xaxis4k, synData_y1_base);
                    y1=ppval(pp1,synData_x_base);
                    y2=ppval(pp2,synData_x_base);
                    
                    MTRasym=y1(end:-1:1)-y1;
                    xbin=zeros(1,4000);
                    xbin(startIdx:stopIdx)=1;
                    asym_integral=(trapz(xbin,MTRasym)); 
                    
%                     %calculate asymSpectrum
%                     asym_Base=abs(y2-y1); 
%                     asym_integral=trapz(asym_Base(startIdx:stopIdx));
                
                    if stopIdx > 4000
                    disp('exceed contrast');
                    stopIdx = max(stopPositions_Contrast);
                    Base_SplineIntegralResults(k,l)=0;
                    else
                    Base_SplineIntegralResults(k,l)=asym_integral;
                    end

                
    end      
end
%  end
 toc
 
cd (path);
% save('mask','mask');
Base_SplineIntegralResults(Base_SplineIntegralResults>100)=0;
figure,imagesc(rot90(Base_SplineIntegralResults,3),[0 0.1]);colormap jet;
colorbar;
title(sprintf('CESTBaseMap%s(spline)',CESTdata.header.patientID));
saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string')));
save('CESTBaseMap(spline).mat','Base_SplineIntegralResults');
save('BaseMap+-5ppm:spline.mat','Base_AcqData','Base_SplineCorrectedData','Base_SplineIntegralResults');
% cd ('/Users/clairezhao/brukermatlabscripts/CEST_DataAnalysis/Matlab_Analysis/InfusionGlucoCEST/');
% ImageViewerStudio(SplineIntegralResults);
end
