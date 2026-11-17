%% load data
clear all;close all;clc;
% CEST part
CESTpath='/Users/clairezhao/Documents/data/pv6/CEST/SEK00-72864__CEST_kidney_02__2019-08-21_10.18.43/E021_wz_RC_FS_128x128_ETL24_E21';
CESTdatapath=fullfile(CESTpath,filesep,'pdata',filesep,'1');
CESTdata=MRIDataLoader.LoadFrom2DSeq(CESTdatapath); 
CEST_Offsets=CESTdata.headerHistory{1}.method.MT_Offsets_NoM0;
CESTdata.Permute('read','line','echo');
handles.CEST.CEST=CESTdata;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


%if no S0data
S0=CESTdata.SubVolume('echo',1);
DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(1);


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
localCEST_S0_Vol = S0.vol;
%handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;
handles.CEST.DataSlopeS0 =DataSlope_S0;
localSizeRead = handles.CEST.CEST.Size('read');
localSizeLine = handles.CEST.CEST.Size('line');

%% pixel-wised analysis with masking
mag_temp=CESTdata(:,:,1);
mask = mag_temp >= median(median(mag_temp)).*20;% approximately thresholding due to experience


 AcqData = zeros(size(localCEST_MT_Offsets_NoM0,2),2,size(mask,1),size(mask,2));
 TLfitData= zeros(2,1000,size(mask,1),size(mask,2));
 SplineData= zeros(2,1000,size(mask,1),size(mask,2));% splined x value,splined y value
 SplineIntegralResults=zeros(size(mask,1),size(mask,2));%integral result
 TLfitResults=zeros(size(mask,1),size(mask,2),1,3);
tic
 for k = 1:localSizeRead
    fprintf('Iteration_Read %3.3i/%3.3i\n', k, localSizeRead); 
    for l = 1:localSizeLine
%             fprintf('Iteration_Line %3.3i/%3.3i\n', l, localSizeLine); 
        if mask(k,l)== 0
            handles.CEST.ContrastMatrixCorrected(k,l) = 0;
%             handles.CEST.RMSEMap(k,l) = 0;
%             handles.CEST.SSEMap(k,l) = 0;
%            handles.CEST.residual(k,l) = 0;
            AcqData(:,:,k,l)=0;
            TLfitData(:,:,k,l)=0;
            TLfitResults(k,l,:,:)=0;
            SplineData(:,:,k,l)=0;
            SplineIntegralResults(k,l)=0;
            continue;
        else
            x = localCEST_MT_Offsets_NoM0;
%             y = ((localCEST_Vol(k,l,:) .* handles.CEST.DataSlope) ./ (localCEST_S0_Vol(k,l,1) .* handles.CEST.DataSlopeS0)).'; 
            y = squeeze((localCEST_Vol(k,l,:))  ./ (localCEST_S0_Vol(k,l,1)))'; 

            [xSort, xIndex] = sort(x);
%                 CEST = [(x(xIndex(:)))*0.0025; y(xIndex(:))]; %ppm
                  CEST_raw = [(x(xIndex(:))); y(xIndex(:))];   % in Hz             
          
                %downsample point for testing
%                 x= downsample(x,5);
%                 y= downsample(y,5);
% -----------------------------direct take minimum ------------------------------------           
                  CEST_corrected=zeros(size(CEST_raw));
%                 CEST_corrected(2,:) = CEST(2,:);
%                   xx = -4:0.01:4;
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
%                 CEST_corrected(1,:) = CEST_raw(1,:) - handles.WASSR.HzOffsetMatrix(k,l);
%                 CEST_corrected(2,:) = CEST_raw(2,:);
%                 
% %               %recheck residual after WASSR correction
%                 [min_val,min_idx]=min(CEST_corrected(2,:));
%                 w0_FromMin=CEST_corrected(1,min_idx);
%                 printf('Detected remaining minimum shift of CEST: %f Hz',w0_FromMin);
% %               %remove residual
%                 CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
%                 handles.CEST.residual(k,l) = w0_FromMin;
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
                AcqData(:,1,k,l)= CEST_corrected(1,:);
                AcqData(:,2,k,l)= CEST_corrected(2,:);
%                 
                
%------------------------------------Multiple Lorentzian-----------------------------------------------------------------
                                 
                %for power calibration series only
%                 [params,xFit,yFit,lorentzInt] = TripleLorentzFit(CEST_corrected(1,1:end-1),CEST_corrected(2,1:end-1));

%               %TripleLorentzFit  
                [params,xFit,yFit,lorentzInt] =CEST.TripleLorentzFit(CEST_corrected(1,:),CEST_corrected(2,:));
                TLfitResults(k,l,:,:)= [lorentzInt(1);lorentzInt(2);lorentzInt(3)];
                TLfitData(:,:,k,l)= [xFit;yFit];
                handles.CEST.ContrastMatrixCorrected(k,l,1,2)= TLfitResults(k,l,1,2)+TLfitResults(k,l,1,3);
                
%-----------------------------------------SPLINE------------------------------------------------------------------------- 
                pp = spline(CEST_corrected(1,:)', CEST_corrected(2,:)');
                Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),1000);
                Spline_Y=ppval(pp,Spline_X);
                SplineData(:,:,k,l)=[Spline_X; ppval(pp,Spline_X)];
                
                start_Contrast = 400;
                startPositions_Contrast = find(Spline_X<start_Contrast);
                startIdx = max(startPositions_Contrast)+1;
%                 
                stop_Contrast = 1200;
                stopPositions_Contrast = find(Spline_X<stop_Contrast);
                stopIdx = max(stopPositions_Contrast)+1;
                
                if stopIdx > 1000
                    disp('exceed contrast');
                    SplineIntegralResults(k,l)=0;
                    
                else
                    
                array2Trapz_x=Spline_X(startIdx:stopIdx);
                array2Trapz_y=1-(Spline_Y(startIdx:stopIdx));%flip upward for integration
                
                SplineIntegralResults(k,l)=trapz(array2Trapz_x,array2Trapz_y);
                end
        end
            
    end
 end
 toc
%% pixel-wise map
 
 
 %% click pop out figure
 h = figure,imagesc(squeeze(SplineIntegralResults(:,:)'),[0 250]),colormap jet;
colormap jet;title('glucoCESTasymSpectrum');
cb = colorbar;
title(cb, 'AUC');

i=0;
while i==0
clear ii jj
[jj,ii]=getpts(h);

ii=int8(ii);
jj=int8(jj);

%generate symmetric datapoint
synData_x=SplineData(1,:,ii,jj);
synData_y=SplineData(2,:,ii,jj);

[syn_min_value,syn_min_idx]=min(synData_y);
synData_y1= synData_y;
FillNr=abs(syn_min_idx-(length(synData_y')./2)).*2;
a=syn_min_idx-(length(synData_y))./2;
if a> 0
synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx-FillNr));
else 
    synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx+FillNr));
end

%calculate asymSpectrum
asym=abs(synData_y-synData_y1);

%plot for single spectrum
figure,plot(AcqData(:,1,ii,jj),AcqData(:,2,ii,jj),'x');
hold on
plot(SplineData(1,:,ii,jj),SplineData(2,:,ii,jj),'b');
axis([-1500 1500 0 1.1]);
plot(synData_x,asym,'k');

end