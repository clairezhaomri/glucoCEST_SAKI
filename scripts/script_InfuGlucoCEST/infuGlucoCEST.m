clear all; close all,clc;
BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1fat/';
cd(BasePath);

% %load highres data
% T2hiresPath=dir('*T2RARE_highres_axialMS*');
% T2hiresDataPath=fullfile(BasePath,T2hiresPath.name,filesep,'pdata',filesep,'1');
% Hiresdata=MRIDataLoader.LoadFrom2DSeq(T2hiresDataPath); 
% Hiresdata.Permute('read','line','echo');


% 
% % CEST InfusionPart(with name on protocol)
% CESTInfuPath=dir('*CEST_infusion*');
% CESTInfuDatapath=fullfile(BasePath,CESTInfuPath(1).name,filesep,'pdata',filesep,'1');
% CESTInfuData=MRIDataLoader.LoadFrom2DSeq(CESTInfuDatapath); 
% CEST_InfuOffsets=CESTInfuData.headerHistory{1}.method.MT_Offsets_NoM0;
% CESTInfuData.Permute('read','line','echo');
% handles.CEST.InfuCEST=CESTInfuData;
% handles.CEST.InfuDataSlope = CESTInfuData.headerHistory{1}.visupars.VisuCoreDataSlope;
cd(BasePath);
load('WASSRMap.mat');
handles.WASSR.WASSRMap=WASSRMap; 


CESTInfuPath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/11/';
CESTInfuDatapath=fullfile(CESTInfuPath,filesep,'pdata',filesep,'1');
CESTInfuData=MRIDataLoader.LoadFrom2DSeq(CESTInfuDatapath); 
CEST_InfuOffsets=CESTInfuData.headerHistory{1}.method.MT_Offsets_NoM0;
CESTInfuData.Permute('read','line','echo');
handles.CEST.InfuCEST=CESTInfuData;
handles.CEST.InfuDataSlope = CESTInfuData.headerHistory{1}.visupars.VisuCoreDataSlope;




%% check M0 location
%------------CEST_infuData---------------------------
if strcmp(CESTInfuData.headerHistory{1}.method.M0Option,'M0_Last')
    InfuS0data=CESTInfuData.vol(:,:,CESTInfuData.Size('echo'));
    CESTInfuData.vol(:,:,CESTInfuData.Size('echo'))=[];
    InfuDataSlope_S0=CESTInfuData.headerHistory{1}.visupars.VisuCoreDataSlope(CESTInfuData.Size('echo'));
    handles.CEST.InfuDataSlope(CESTInfuData.Size('echo')) = [];
end

if strcmp(CESTInfuData.headerHistory{1}.method.M0Option,'M0_First')
    InfuS0data=CESTInfuData.vol(:,:,1);
    CESTInfuData.vol(:,:,1)=[];
    InfuDataSlope_S0=CESTInfuData.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    handles.CEST.InfuDataSlope(1) = [];
end

InfuS0data = MRIDataBasic(InfuS0data, 'read,line,slice,echo');

%% makedir local
localInfuCEST_MT_Offsets_NoM0 = handles.CEST.InfuCEST.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localInfuCEST_Vol = handles.CEST.InfuCEST.vol;
localInfuCEST_S0_Vol = InfuS0data;
handles.CEST.InfuDataSlope_S0=InfuDataSlope_S0;
localSizeRead = handles.CEST.InfuCEST.Size('read');
localSizeLine = handles.CEST.InfuCEST.Size('line');
%% pixelwise infusion map

load('mask.mat');

 Infu_AcqData = zeros(size(localInfuCEST_MT_Offsets_NoM0,2),2,size(mask,1),size(mask,2));
 %TLfitData= zeros(2,1000,size(mask,1),size(mask,2));
 Infu_SplineCorrectedData= zeros(size(mask,1),size(mask,2),2,100);% splined x value,splined y value
 Infu_SplineIntegralResults=zeros(size(mask,1),size(mask,2));%integral result
 %TLfitResults=zeros(size(mask,1),size(mask,2),1,3);
 
tic
 for k = 1:localSizeRead
    fprintf('Iteration_Read %3.3i/%3.3i\n', k, localSizeRead); 
    for l = 1:localSizeLine
%             fprintf('Iteration_Line %3.3i/%3.3i\n', l, localSizeLine); 
        if mask(k,l)== 0
%             handles.CEST.ContrastMatrixCorrected(k,l) = 0;
%             handles.CEST.RMSEMap(k,l) = 0;
%             handles.CEST.SSEMap(k,l) = 0;
%            handles.CEST.residual(k,l) = 0;
            Infu_AcqData(:,:,k,l)=0;
            %TLfitData(:,:,k,l)=0;
            %TLfitResults(k,l,:,:)=0;
            Infu_SplineCorrectedData(k,l,:,:)=0;
            Infu_SplineIntegralResults(k,l)=0;
            continue;
        else
            x = localInfuCEST_MT_Offsets_NoM0;
%             y = ((localInfuCEST_Vol(k,l,:) .* handles.CEST.InfuDataSlope') ./ (localInfuCEST_S0_Vol(k,l,1) .* handles.CEST.InfuDataSlope_S0)); 
            y = squeeze((localInfuCEST_Vol(k,l,:))  ./ (localInfuCEST_S0_Vol(k,l,1)))'; 

            [xSort, xIndex] = sort(x);
%                 CEST = [(x(xIndex(:)))*0.0025; y(xIndex(:))]; %ppm
                  CEST_raw = [(x(xIndex(:))); y(xIndex(:))];   % in Hz             
          
                %downsample point for testing
%                 x= downsample(x,5);
%                 y= downsample(y,5);
% -----------------------------direct take minimum ------------------------------------           
%                   CEST_corrected=zeros(size(CEST_raw));
% %                 CEST_corrected(2,:) = CEST(2,:);
% %                   xx = -4:0.01:4;
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
                Infu_AcqData(:,1,k,l)= CEST_corrected(1,:);
                Infu_AcqData(:,2,k,l)= CEST_corrected(2,:);
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
                Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),100);
                Spline_Y=ppval(pp,Spline_X);
                [min_val,min_idx]=min(Spline_Y);
                w0_FromSplineMin=Spline_X(min_idx);
                Spline_Xcorrected=Spline_X - w0_FromSplineMin;
                Infu_SplineCorrectedData(k,l,:,:)=[Spline_Xcorrected; ppval(pp,Spline_Xcorrected)];
                %update AcqData x axis
                Infu_AcqData(:,1,k,l)= CEST_corrected(1,:)-(w0_FromSplineMin./100);
                
                
                start_Contrast = 400;
                startPositions_Contrast = find(Spline_X<start_Contrast);
                startIdx = max(startPositions_Contrast)+1;
%                 
                stop_Contrast = 1200;
                stopPositions_Contrast = find(Spline_X<stop_Contrast);
                stopIdx = max(stopPositions_Contrast)+1;
                
                if stopIdx > 100
                    disp('exceed contrast');
                    Infu_SplineIntegralResults(k,l)=0;
                    
                else
                    
                array2Trapz_x=Spline_X(startIdx:stopIdx);
                array2Trapz_y=1-(Spline_Y(startIdx:stopIdx));%flip upward for integration
                
                Infu_SplineIntegralResults(k,l)=trapz(array2Trapz_x,array2Trapz_y);
                end
        end
            
    end
 end
 toc
 
cd (BasePath);
% uisave({'Infu_AcqData','Infu_SplineCorrectedData','Infu_SplineIntegralResults'},'InfuMap');
save('InfuMap100.mat','Infu_AcqData','Infu_SplineCorrectedData','Infu_SplineIntegralResults');

%% %% quality assurance pop up figure
% % h = figure,imagesc(rot90(Infu_SplineIntegralResults(:,:),3),[0 350]),colormap jet;
% h = figure,imagesc(Infu_SplineIntegralResults(:,:),[0 350]),colormap jet; % rotate 90 of every animal dataset, was not sure how to turn it yet.
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
% %generate symmetric datapoint(only works in phantom!!)
% synData_x=Infu_SplineCorrectedData(ii,jj,1,:);
% synData_y=Infu_SplineCorrectedData(ii,jj,2,:);
% 
% [syn_min_value,syn_min_idx]=min(synData_y);
% synData_y1= synData_y;
% FillNr=abs(syn_min_idx-(length(synData_y)./2)).*2;
% a=syn_min_idx-(length(synData_y))./2;
% if a> 0
% synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx-FillNr));
% else 
%     synData_y1(syn_min_idx+1:end)=fliplr(synData_y(1:syn_min_idx+FillNr));
% end
% 
% %calculate asymSpectrum
% asym=abs(synData_y-synData_y1);
% 
% %plot for single spectrum
% figure,plot(Infu_AcqData(:,1,ii,jj),Infu_AcqData(:,2,ii,jj),'x');
% hold on
% plot(squeeze(Infu_SplineCorrectedData(ii,jj,1,:)),squeeze(Infu_SplineCorrectedData(ii,jj,2,:)),'b');
% axis([-1800 1800 0 1.1]);
% % % plot(synData_x(startIdx:stopIdx),asym(startIdx:stopIdx),'k');
% % plot(squeeze(synData_x(startIdx:stopIdx)),squeeze(asym(startIdx:stopIdx)),'k');
% 
% end