function [Data,GlobalM0Map,T1Map_lsq,T2Map] = T1T2MapFitting(DataPath)

%% T1Fitting for Bruker Data

DataPath

%figure

% load Bruker 2dseq data (full path)
Data=MRIDataLoader.LoadFrom2DSeq([DataPath filesep 'pdata' filesep '1' filesep '2dseq']);

% check structure: read,phase,slice,global/selective IR,16 different TI times   
size(Data.vol)

[sval,sindx]=sort(Data.headerHistory{1}.acqp.ACQ_recov_time);

Data.headerHistory{1}.acqp.ACQ_recov_time=Data.headerHistory{1}.acqp.ACQ_recov_time(sindx);
Data.vol=Data.vol(:,:,:,[sindx]);

if (false) % roi entry.... 
myHandle=ImageViewer(Data);
waitfor(myHandle.figure);

roiData=Data.Copy();

roiData.Linearize('read,line');
roiData.Permute('linearized');
%data.Mask=Data.SubVolume('echo',2);

for iRoi = 1:size(myHandle.ROIs,1)
    data.rois{iRoi}.roimask=myHandle.ROIs(iRoi).volume.vol;
    data.rois{iRoi}.roiname=myHandle.ROIs(iRoi).name;
    data.rois{iRoi}.roisize=length(find(data.rois{iRoi}.roimask(:)));
%    data.Mask.vol=data.Mask.vol+data.rois{iRoi}.roimask.*iRoi;
    % store the data values in each ROIs  (mean, std)
    data.roiVals{iRoi}.mean=mean(roiData(data.rois{iRoi}.roimask,:,:,:,:,:,:),'omitnan');
    data.roiVals{iRoi}.std=std(roiData(data.rois{iRoi}.roimask,:,:,:,:,:,:),'omitnan');
end
end

% Martins Democode
% fix Dataume dimensions

Data.Permute('read,line,slice,echo,channel');

% vdim = VolumeDimensions;
% if ndims(Data.vol) == 5    
%     Data.dimType.Create([vdim.read vdim.line vdim.slice vdim.ida vdim.echo]);
% else
%     Data.dimType.Create([vdim.read vdim.line vdim.ida vdim.echo]);
% end
% 
% Data.Permute([vdim.read vdim.line vdim.slice vdim.ida vdim.echo],'Squeeze',false);
% 
% OrigDims=size(Data.vol);

% % smooth with savitzky golay
%DataSmooth = Data.Copy;
% DataSmooth.Smooth(, 7);
% 
% size(DataSmooth.vol)

% 
% % smooth with imfilter
% H = fspecial('gaussian',[10 10], 2);
% % DataBlur = Data.Copy();
% DataSmooth = imfilter(Data,H,'replicate');

% Data=reshape(DataSmooth.vol,OrigDims);
% DataSmooth.SaveToFile([DataPath filesep  'DataSmooth.nii']);

% if ndims(Data) == 5    
% disp 'Data is a 5D array  -- x,y,z,2,TIs -- good!'
% 
% elseif  ndims(Data) == 4    
%     Data=reshape(Data,size(Data,1),size(Data,2),1,size(Data,3),size(Data,4));
% else
%     disp 'Data should be a 5D array  but ISN´T -- x,y,z,2,TIs '
%     size(Data)
%     
% end


%% Filter Data
% funktioniert nur mit 2D Daten
%Data=medfilt2(Data,[3 3]);

%% Slice loop 


TI=transpose(Data.headerHistory{1}.acqp.ACQ_recov_time);
% TI=Data.headerHistory{1}.acqp.ACQ_recov_time;
EffectiveTE=Data.headerHistory{1}.acqp.ACQ_echo_time;
EffectiveTE=EffectiveTE';

if size(TI,1) ~= size(Data,5)
    error 'Error: Size mismatch TI times'
end
%% for motion data drop before looping
%  Data(:,:,:,:,4)=[];
% % % % % Data(:,:,:,:,1)=[];
% % % % % Data(:,:,:,:,1)=[];
% % % % % % Data(:,:,:,:,3)=[];
%  TI(4)=[];
% % % TI(1)=[];% delete 1st and 2nd motion echoes if necessary
% 
% % % TI(3)=[];
% % % TI(9)=[];
%% z,x,y loop for fitting 
%prepare Data with masking.
cd (DataPath);

if isfile('mask.mat')
   load('mask.mat');
else
    ssImg=Data(:,:,1,1,3);% single slice Image
    I=figure,imagesc(ssImg);
    mask=roipoly;
    waitfor(I);
    mask=repmat(mask(:,:),[1 1 size(Data,3) size(Data,4) size(Data,5)]);
end

Data=Data.* mask;

save('mask','mask');
save('Data','Data');
cd('/Users/clairezhao/Documents/MATLAB/CEST/script_T1Map/');

%drop data before looping
%  TI(2)=[];
%  TI(1)=[];
%  Data(:,:,:,:,2)=[];
%  Data(:,:,:,:,1)=[];


T1Map_lsq=zeros(size(Data,1),size(Data,2), size(Data,3));
% GlobalMap=zeros(size(Data,1),size(Data,2), size(Data,3));
% MSEMap=zeros(size(Data,1),size(Data,2), size(Data,3));
fitResults=zeros(3,1,size(Data,1),size(Data,2), size(Data,3));
AcqData=zeros(length(TI),1,size(Data,1),size(Data,2));
% T2Map=zeros(size(Data,1),size(Data,2), size(Data,3));
% T2fitResults=zeros(3,1,size(Data,1),size(Data,2), size(Data,3));
% T2AcqData=zeros(length(EffectiveTE),1,size(Data,1),size(Data,2));

T1ResidualMap=zeros(size(Data,1),size(Data,2), size(Data,3));
% Data(:,:,:,:,1:2)=Data(:,:,:,:,1:2).*(-1.0);

tic
for sl=1:size(Data,3)
    for xi=1:size(Data,1)
        for yi=1:size(Data,2)
            % skip if outside mask
            if mask(xi, yi, sl,1,1)== 0
                T1Map_lsq(xi, yi, sl) = 0;  
%                 GlobalMap(xi, yi, sl) = 0;                
                T1ResidualMap(xi, yi, sl) = 0;
                fitResults(3,1,xi,yi,sl)=0;
% %                 MSEMap(xi,yi,sl)=0;

%                     T2Map(xi, yi, sl) = 0;  
%                     T2fitResults(:,1,xi,yi,sl)=0;
%                    T2AcqData(:,1,xi,yi,sl)=0;
                continue;
            else
 %--------------- fit T1 Vector----------------------
                FitVector=(squeeze(Data(xi,yi,sl,2,:)));
%                   FitVector_smooth=smooth(FitVector);
 %               subplot(2,1,1); 
%                plot(TI,FitVector,'-o'); 
                try
                    options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
                    lb=[0 0 0];
                    ub=[];% 
                    A0=[50,1000,1];
%                     fun=@(A,X)A(1)*(1-2*exp(-X./A(2)))+A(3);
                    fun=@(A,X)A(1).*(1-(2-A(3)).*exp(-X./A(2)));%brukerfunction
                    [Params,residual] = lsqcurvefit(fun, A0, TI, FitVector,lb,ub,options);
%                     Params = nlinfit(TI, FitVector_smooth, @T1Map_Fit, [6 1200 9],'ErrorModel','combined');
%                     GlobalM0Map(xi, yi, sl) = Params(1);%signal level
                    T1Map_lsq(xi, yi, sl) = Params(2);%T1 Relaxation Time
                    fitResults(:,1,xi,yi,sl)=[Params(1) Params(2) Params(3)];
                    AcqData(:,1,xi,yi,sl)=FitVector;
                    T1ResidualMap(xi,yi,sl)=residual;
%                      MSEMap(xi, yi, sl) = MSE;
                    % Debugging Plots
 %                   figure,plot(TI,FitVector,'o');
%                    hold on
%                    plot(linspace(TI(1),TI(end),100),T1Map_Fit(Params,linspace(TI(1),TI(end),100)),'-'); 
%                     hold off
                 catch Me
                    TI';
                    FitVector';
 %                   plot(TI,FitVector)
 %                   pause(0.5)
%                    GlobalM0Map(xi, yi, sl) = 0;
                    T1Map_lsq(xi, yi, sl) = 0;
                    fitResults(:,1,xi,yi,sl)=0;
                    AcqData(:,1,xi,yi,sl)=FitVector;
                    T1ResidualMap(xi,yi,sl)=0;
%                     save('T1Map','T1Map');
                end

%--------------------- fit T2 Vector------------------
%                 FitVector=squeeze(Data(xi,yi,sl,:,1)); 
%  %               subplot(2,1,2); plot(TI,FitVector,'-o');
%                 try
%                     Params = nlinfit(EffectiveTE, FitVector, @T2Map_Fit, [65 29 9]);
%                     T2fitResults(:,1,xi,yi,sl)=[Params(1) Params(2) Params(3)];
%                     T2AcqData(:,1,xi,yi,sl)=FitVector;
% %                     GlobalM0Map(xi, yi, sl) = Params(1);
%                     T2Map(xi, yi, sl) = Params(2);
%                     % Debugging Plots
% %                    hold on
% %                    plot(linspace(TI(1),TI(end),100),T1_Fit_FAIR(Params,linspace(TI(1),TI(end),100)),'-'); 
% %                    hold off
%                 catch Me
%                     EffectiveTE'
%                     FitVector'
% %                    plot(TI,FitVector)
% %                    pause(0.5)
% %                     SelectiveM0Map(xi, yi, sl) = 0;
%                     T2Map(xi, yi, sl) = 0;
%                     T2fitResults(:,1,xi,yi,sl)=0;
%                     T2AcqData(:,1,xi,yi,sl)=0;
%                 end
% %                pause(0.1 ); % show updated debugging plots 
% %              else
% % %                GlobalM0Map(xi, yi, sl) = 0;
% %                T2Map(xi, yi, sl) = 0;
%-----------------------------------------------------
            

            end

       end
    
   end
end

toc
   cd(DataPath);
   save('mask','mask');
%    ImageViewerStudio(T1Map);
% save('T2AcqDataAndFit.mat','T2Map', 'T2fitResults','T2AcqData','EffectiveTE');

save('T1Map_lsq','T1Map_lsq');
save('T1ResidualMap','T1ResidualMap');
save('AcqAndLsqFit.mat','AcqData','fitResults','TI');
% MRIDataSaver.SaveToNIFTI('T1Map.nii',T1Map);
% save('T2Map','T2Map');
% save('GlobalM0Map','GlobalM0Map');
MRIDataSaver.SaveToNIFTI('T1ResidualMap.nii',T1ResidualMap);
MRIDataSaver.SaveToNIFTI('T1Map_lsq.nii',T1Map_lsq);
% MRIDataSaver.SaveToNIFTI('T2Map.nii',T2Map);
% MRIDataSaver.SaveToNIFTI('GlobalM0Map.nii',GlobalM0Map);
% % 
% % save('T2Map','T2Map');
% % save('GlobalM0Map','GlobalM0Map');



% %debug plot
% h = figure,imagesc(squeeze(T1Map_lsq(:,:,1)),[0 1000]),colormap jet;
% colormap jet;title('T1Map:lsq');
% cb = colorbar;
% title(cb, 'ms');
% 
% i=0;
% while i==0
% clear ii jj
% [jj,ii]=getpts(h);
% 
% ii=int8(ii);
% jj=int8(jj);
% 
% %plot for single spectrum
% figure,plot(TI,AcqData(:,1,ii,jj),'o');
% hold on
% fitResultY=fun(fitResults(:,1,ii,jj),linspace(TI(1),TI(end),100));
% plot(linspace(TI(1),TI(end),100),fitResultY,'-');
% title(sprintf('wt#71983_1_0_0_4Base @(%d,%d)',ii,jj));
% end

cd('/Users/clairezhao/Documents/MATLAB/CEST/script_T1Map');
end

