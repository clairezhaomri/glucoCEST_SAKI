function [Data,WASSRMap,WASSR_AcqData, FWHMMap]= WASSR(DataPath)
%% load WASSR data
DataPath

Data=MRIDataLoader.LoadFrom2DSeq([DataPath filesep 'pdata' filesep '1' filesep '2dseq']);



% clear all;close all;clc;
% Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/';
% WASSRpath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200107_190957_MIBIO_00465_1_1/10/';
% WASSRdatapath=fullfile(WASSRpath,filesep,'pdata',filesep,'1');
% Data=MRIDataLoader.LoadFrom2DSeq(WASSRdatapath); 
WASSR_Offsets=Data.headerHistory{1}.method.MT_Offsets_NoM0;
Data.Permute('echo','read','line');
handles.WASSR.WASSR=Data;
handles.WASSR.DataSlope = Data.headerHistory{1}.visupars.VisuCoreDataSlope;



% % external S0
% WS0path='/Users/clairezhao/Documents/data/MRI/pv6/CEST/vh170107.Fq1/24WASSRS0/';
% WS0datapath=fullfile(WS0path,filesep,'pdata',filesep,'1');
% WS0data=MRIDataLoader.LoadFrom2DSeq(WS0datapath); 
% WS0data.Permute('read','line','echo');
% DataSlope_S0=WS0data.headerHistory{1}.visupars.VisuCoreDataSlope(1);

%% check M0 location

if strcmp(Data.headerHistory{1}.method.M0Option,'M0_Last')
    WS0data=Data.vol(Data.Size('echo'),:,:);
    Data.vol(Data.Size('echo'),:,:)=[];
    DataSlope_S0=Data.headerHistory{1}.visupars.VisuCoreDataSlope(Data.Size('echo'));
    handles.WASSR.DataSlope(Data.Size('echo')) = [];
end

if strcmp(Data.headerHistory{1}.method.M0Option,'M0_First')
    WS0data=Data.vol(1,:,:);
    Data.vol(1,:,:)=[];
    DataSlope_S0=Data.headerHistory{1}.visupars.VisuCoreDataSlope(1);
    handles.WASSR.DataSlope(1) = [];
end

%% make directory local


[xSortW, xIndexW] = sort(WASSR_Offsets);
Wlocscaleddata=handles.WASSR.WASSR.vol ./ WS0data;
Wlocscaleddata=Wlocscaleddata(xIndexW,:,:);

localWASSR_Offsets_NoM0 = handles.WASSR.WASSR.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localWASSR_Vol = handles.WASSR.WASSR.vol;
localWASSR_S0_Vol = WS0data;

localSizeRead = handles.WASSR.WASSR.Size('read');
localSizeLine = handles.WASSR.WASSR.Size('line');%%



  %% Roi-based analysis and plot
% cd (DataPath);
% if isfile('ROI.mat')
%    load('ROI.mat');
% else
%     I=figure,imagesc(squeeze(Data(1,:,:)));colormap jet;
%     ROI=roipoly;
%     waitfor(I);
%     save('ROI','ROI');
% end
%  
%  
% RoiMeanW=zeros(handles.WASSR.WASSR.headerHistory{1}.method.Number_Offsets,1);
% RoiStdW=zeros(handles.WASSR.WASSR.headerHistory{1}.method.Number_Offsets,1);
% 
% for echoi = 1:handles.WASSR.WASSR.headerHistory{1}.method.Number_Offsets
%    singleecho=Wlocscaleddata(echoi,:,:);
%    RoiMeanW(echoi)=mean2(singleecho(ROI));
%    RoiStdW(echoi)=std2(singleecho(ROI));
% end
% 
% figure;plot (xSortW,RoiMeanW,'-*');
% hold on;
% % plot (xSortW,RoiStdW.*10);
% % hold off;
% 
% pause(1);

%% 
cd (DataPath);
if isfile('mask.mat')
   load('mask.mat');
else
    mag_temp=squeeze(Data(1,:,:));
    I=figure,imagesc(mag_temp(:,:));
    mask=roipoly;
    waitfor(I);
    save('mask','mask');
end

WASSRMap=zeros(size(mask,1),size(mask,2));
WASSR_AcqData = zeros(size(localWASSR_Offsets_NoM0,2),2,size(mask,1),size(mask,2));
FWHMMap=zeros(size(mask,1),size(mask,2));
 
tic
for k=1:localSizeRead
    fprintf('Iteration_Line %3.3i/%3.3i\n', k, localSizeRead); 
    for l=1:localSizeLine
        if mask(k,l)== 0
             WASSRMap(k,l) = 0;
             WASSR_AcqData(:,:,k,l)=0;
             FWHMMap(k,l)=0;
        continue;
    else
        %code for calculating the spectrum
                 x = localWASSR_Offsets_NoM0;
                 y = ((localWASSR_Vol(:,k,l) .* handles.WASSR.DataSlope.') ./ (localWASSR_S0_Vol(1,k,l) .* DataSlope_S0)).'; 

                 [sval,sindx]=sort(x);
                 x_sorted=x(sindx(:));
                 y_sorted=y(sindx(:));
               
%% Fit: 'Lorentzian Fit'.
                [xData, yData] = prepareCurveData( x_sorted, y_sorted );

                % Set up fittype and options.
                ft = fittype( '(d-((a .*(1./(pi.* c))) ./ (1+ (b-x).^2 .* (1./(pi.*c)).^2)))', 'independent', 'x', 'dependent', 'y' );
                %a: d:as a inverse of 1,
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.StartPoint = [0.400271854321867 0.189264381642447 0.928792973953908 0.957166948242946];

                % Fit model to data.
                [fitresult, gof] = fit( xData, yData, ft, opts );

                % debug Plot fit with data.
%                 figure( 'Name', 'Lorentzian Fit' );
%                 h = plot( fitresult, xData, yData );
%                 legend( h, 'y_sorted vs. x_sorted', 'Lorentzian Fit', 'Location', 'NorthEast' );
%                 % Label axes
%                 xlabel x_sorted
%                 ylabel y_sorted
%                 grid on

                 WASSR_AcqData(:,1,k,l)= x_sorted;
                 WASSR_AcqData(:,2,k,l)= y_sorted;
%                  LfitData(1,:,k,l)= xx;
%                  LfitData(2,:,k,l)= fun(Params,xx)';
%                  LfitResults(k,l,:,:)=Params;
                 WASSRMap(k,l) = fitresult.b;
                 FWHMMap(k,l) = (1./(pi.*(fitresult.c)));

%------------------------------------spline-----------------------------------------------------------
%                     qq = spline(x_sorted(1,:), y_sorted(1,:));
%        
%                     yMinValue = min(ppval(qq,xx));
%                     xPos = find(ppval(qq,xx) == yMinValue);
%                     xValue = xx(xPos);
%                     if size(xValue) == [1 1]
%                         handles.WASSR.OffsetMatrix(k,l) = xValue;
%                     else
%                         disp 'wierd xavlues size'
%                         xValue
%                         size(xValue)
%                     handles.WASSR.OffsetMatrix(k,l) =0;
%                     end
% 
%------------------------------------------------------------------------------------------------------
         end 
     end
end
toc

save('WASSRMap','WASSRMap');
figure,imagesc(WASSRMap,[-300 300]);
colorbar
title(sprintf('WASSRMap:%s',Data.header.patientID));
saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));

save('WASSRResult.mat','WASSR_AcqData','WASSRMap','FWHMMap');


cd('/Users/clairezhao/Documents/MATLAB/CEST/script_InfuGlucoCEST');
end

