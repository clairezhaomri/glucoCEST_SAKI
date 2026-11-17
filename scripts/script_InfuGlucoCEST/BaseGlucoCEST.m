clear all,close all,clc;
%% load data
% BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/ANEST-71982__HemeCEST_06_Heme__2019-10-07_19.38.43/';
% cd(BasePath);

% T2hiresPath=dir('*T2RARE_highres_axialMS*');
% T2hiresDataPath=fullfile(BasePath,T2hiresPath.name,filesep,'pdata',filesep,'1');
% Hiresdata=MRIDataLoader.LoadFrom2DSeq(T2hiresDataPath); 
% Hiresdata.Permute('read','line','echo');

% CEST part
% CESTpath=dir('*CEST_baseline*');
% CESTdatapath=fullfile(BasePath,CESTpath.name,filesep,'pdata',filesep,'1');
% CESTdata=MRIDataLoader.LoadFrom2DSeq(CESTdatapath); 
% CEST_Offsets=CESTdata.headerHistory{1}.method.MT_Offsets_NoM0;
% CESTdata.Permute('read','line','echo');
% handles.CEST.CEST=CESTdata;
% handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


%bruker path in number
BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/';
CESTpath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/11/';
CESTdatapath=fullfile(CESTpath,filesep,'pdata',filesep,'1');
CESTdata=MRIDataLoader.LoadFrom2DSeq(CESTdatapath); 
CEST_Offsets=CESTdata.headerHistory{1}.method.MT_Offsets_NoM0;
CESTdata.Permute('read','line','echo');
handles.CEST.CEST=CESTdata;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


cd(BasePath);
load('WASSRMap.mat');
handles.WASSR.WASSRMap=WASSRMap; 
% 
% 
% % external S0
% S0path='/Users/clairezhao/Documents/data/MRI/pv6/CEST/vh170107.Fq1/27S0/';
% S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
% S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath); 
% S0data.Permute('read','line','echo');
% DataSlope_S0=S0data.headerHistory{1}.visupars.VisuCoreDataSlope(1);


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

% S0data = MRIDataBasic(S0data, 'read,line,slice,echo');

%% make directory local
localCEST_MT_Offsets_NoM0 = handles.CEST.CEST.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localCEST_Vol = handles.CEST.CEST.vol;
localCEST_S0_Vol = S0data;
handles.CEST.DataSlopeS0 =DataSlope_S0;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


localSizeRead = handles.CEST.CEST.Size('read');
localSizeLine = handles.CEST.CEST.Size('line');

%%
cd (BasePath);
if isfile('mask.mat')
   load('mask.mat');
else
    mag_temp=CESTdata(:,:,1);
    I=figure,imagesc(mag_temp(:,:));
    mask=roipoly;
    waitfor(I);
    save('mask','mask');
end

 Base_AcqData = zeros(size(localCEST_MT_Offsets_NoM0,2),2,size(mask,1),size(mask,2));
 %TLfitData= zeros(2,1000,size(mask,1),size(mask,2));
 Base_SplineCorrectedData= zeros(size(mask,1),size(mask,2),2,100);% splined x value,splined y value
 Base_SplineIntegralResults=zeros(size(mask,1),size(mask,2));%integral result
 %TLfitResults=zeros(size(mask,1),size(mask,2),1,3);
 
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
            Base_AcqData(:,:,k,l)=0;
            %TLfitData(:,:,k,l)=0;
            %TLfitResults(k,l,:,:)=0;
            Base_SplineCorrectedData(k,l,:,:)=0;
            Base_SplineIntegralResults(k,l)=0;
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
% %                 CEST_corrected(2,:) = CEST(2,:);
% %                   xx = -4:0.01:4;
%                   CEST_corrected=zeros(size(CEST_raw));
%                   [min_val,min_idx]=min(CEST_raw(2,:));
%                   w0_FromMin=CEST_raw(1,min_idx);
% %                 printf('Detected remaining minimum shift of CEST: %f ppm',w0_FromMin);
%                 %remove residual
%                   CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
%                   CEST_corrected(2,:) = CEST_raw(2,:);
%                   
% --------------------%for power calibration series ONLY --------------------------- 
%                    CEST_corrected(2,:) = CEST_corrected(2,:)./CEST_corrected(2,end);

%--------------------------------------Use WASSR----------------------------------
%                 CEST_corrected=zeros(size(CEST_raw));
                CEST_corrected(1,:) = CEST_raw(1,:) - handles.WASSR.WASSRMap(k,l);
                CEST_corrected(2,:) = CEST_raw(2,:);
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
                pp = spline(CEST_corrected(1,:)', CEST_corrected(2,:)');
                Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),100);% test effect of 100spline and 1000spline
                Spline_Y=ppval(pp,Spline_X);
                [min_val,min_idx]=min(Spline_Y);
                w0_FromSplineMin=Spline_X(min_idx);
                Spline_Xcorrected=Spline_X - w0_FromSplineMin;
                Base_SplineCorrectedData(k,l,:,:)=[Spline_Xcorrected; ppval(pp,Spline_Xcorrected)];
                %update AcqData x axis
                Base_AcqData(:,1,k,l)= CEST_corrected(1,:)-(w0_FromSplineMin./100);
                
                
                start_Contrast = 400;
                startPositions_Contrast = find(Spline_X<start_Contrast);
                startIdx = max(startPositions_Contrast)+1;
%                 
                stop_Contrast = 1200;
                stopPositions_Contrast = find(Spline_X<stop_Contrast);
                stopIdx = max(stopPositions_Contrast)+1;
                
                if stopIdx > 100
                    disp('exceed contrast');
                    Base_SplineIntegralResults(k,l)=0;
                    
                else
                    
                array2Trapz_x=Spline_X(startIdx:stopIdx);
                array2Trapz_y=1-(Spline_Y(startIdx:stopIdx));%flip upward for integration
                
                Base_SplineIntegralResults(k,l)=trapz(array2Trapz_x,array2Trapz_y);
                end
        end
            
    end
 end
 toc
 
cd (BasePath);
save('mask','mask');
save('BaseMap100.mat','Base_AcqData','Base_SplineCorrectedData','Base_SplineIntegralResults');

% ImageViewerStudio(SplineIntegralResults);

%% quality assurance pop up figure
% % h = figure,imagesc(rot90(Infu_SplineIntegralResults(:,:),3),[0 350]),colormap jet;
% h = figure,imagesc(Base_SplineIntegralResults(:,:),[0 500]),colormap jet; % rotate 90 of every animal dataset, was not sure how to turn it yet.
% 
% colormap jet;title('glucoCESTasymSpectrum');
% cb = colorbar;
% title(cb, 'AUC');
% 
% i=0;
% while i==0
% clear ii jj
% [jj,ii]=getpts(h);
% 
% ii=int8(ii);
% jj=int8(jj);
% 
% % %generate symmetric datapoint(only works in phantom!!)
% % synData_x=Base_SplineCorrectedData(ii,jj,1,:);
% % synData_y=Base_SplineCorrectedData(ii,jj,2,:);
% % 
% % [syn_min_value,syn_min_idx]=min(synData_y);
% % synData_y1= synData_y;
% % FillNr=abs(syn_min_idx-(length(synData_y)./2)).*2;
% % a=syn_min_idx-(length(synData_y))./2;
% % if a> 0
% % synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx-FillNr));
% % else 
% %     synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx+FillNr));
% % end
% % 
% % %calculate asymSpectrum
% % asym=abs(synData_y-synData_y1);
% 
% %plot for single spectrum
% figure,plot(Base_AcqData(:,1,ii,jj),Base_AcqData(:,2,ii,jj),'x');
% hold on
% plot(squeeze(Base_SplineCorrectedData(ii,jj,1,:)),squeeze(Base_SplineCorrectedData(ii,jj,2,:)),'b');
% axis([-1200 1200 0 1.1]);
% % % plot(synData_x(startIdx:stopIdx),asym(startIdx:stopIdx),'k');
% % plot(squeeze(synData_x(startIdx:stopIdx)),squeeze(asym(startIdx:stopIdx)),'k');
% 
% end