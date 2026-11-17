function [AnimalID,SelectFreqAsymResult] = SelectFreqAsymCESTMap(path,SelectPpm)
%Wan-Ting Zhao. Univerisitätsklinikum Jena
%Version 2020 Oct. 13
%This script spline and direct take minimum pixelwisely.


       CESTdata=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']);
       AnimalID = CESTdata.header.patientID;
%        ACQ_scan_name = CESTdata.headerHistory{1}.acqp.ACQ_scan_name;
%        ACQ_time = CESTdata.headerHistory{1}.acqp.ACQ_time;

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


%  Base_AcqData = zeros(size(localCEST_MT_Offsets_NoM0,2),2,128,128);
%  Base_SplineCorrectedData= zeros(128,128,2,4000);% splined x value,splined y value
 SelectFreqAsymResult=zeros(size(CESTdata,1),size(CESTdata,2));%integral result
 
%pixelwise spline and store value  
tic
 for k = 1:localSizeRead
    fprintf('Iteration_Read %3.3i/%3.3i\n', k, localSizeRead); 
    for l = 1:localSizeLine
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
          

% -----------------------------direct take minimum ------------------------------------           
% % 
                  CEST_corrected=zeros(size(CEST_raw));
                  [min_val,min_idx]=min(CEST_raw(2,:));
                  w0_FromMin=CEST_raw(1,min_idx);
%                 printf('Detected remaining minimum shift of CEST: %f ppm',w0_FromMin);
                %remove residual
                  CEST_corrected(1,:) = CEST_raw(1,:) - w0_FromMin;
                  CEST_corrected(2,:) = CEST_raw(2,:);

%-----------------------------------------SPLINE------------------------------------------------------------------------- 
%                 CEST_corrected_X=(CEST_corrected(1,:)');
%                 CEST_corrected_Y=(CEST_corrected(2,:)');

                CEST_corrected_Spline_X=linspace(CEST_corrected(1,1),CEST_corrected(1,end),4000);                
                Spline_Y=spline(CEST_corrected(1,:),CEST_corrected(2,:),CEST_corrected_Spline_X);

                
%                 %debug plot
%                 figure,plot(CEST_corrected_X,CEST_corrected_Y,'b*');
%                 hold on;
%                 plot(CEST_corrected_Spline_X, Spline_Y,'r-');
%                 hold off;

                [min_val,min_idx]=min(Spline_Y);
                w0_FromSplineMin=CEST_corrected_Spline_X(min_idx);
                Spline_Xcorrected=CEST_corrected_Spline_X - w0_FromSplineMin;
                
                 pp1=spline(Spline_Xcorrected,Spline_Y);
                 xaxis=linspace(-2000, 2000, 4000);
                 splinedROIMean4k=ppval(pp1,xaxis);%splineData
                
                SelectFreq=SelectPpm*400;
                glucoCESTSpectumValue_Positive=ppval(pp1,SelectFreq);
                glucoCESTSpectumValue_Negative=ppval(pp1,-SelectFreq);
                
                if SelectPpm==3.5%APT
                GlucoCESTWorkHz4ppm=4*400;
                GlucoCESTWorkHz3ppm=3*400;

                glucoCESTSpectumValue_4ppm=ppval(pp1,GlucoCESTWorkHz4ppm);
                glucoCESTSpectumValue_3ppm=ppval(pp1,GlucoCESTWorkHz3ppm);
                glucoCESTSpectumValue_APT=ppval(pp1,SelectFreq);
                Asym_SelectFreq=(glucoCESTSpectumValue_4ppm+glucoCESTSpectumValue_3ppm)./2-glucoCESTSpectumValue_APT;

                elseif SelectPpm==-3.5
                    GlucoCESTWorkHzneg2ppm=-2*400;
                    GlucoCESTWorkHz5ppm=5*400;

                    glucoCESTSpectumValue_neg2ppm=ppval(pp1,GlucoCESTWorkHzneg2ppm);
                    glucoCESTSpectumValue_5ppm=ppval(pp1,GlucoCESTWorkHz5ppm);
                    glucoCESTSpectumValue_NOE=ppval(pp1,SelectPpm);
                    Asym_SelectFreq=(glucoCESTSpectumValue_neg2ppm+glucoCESTSpectumValue_5ppm)./2-glucoCESTSpectumValue_NOE;

                else
                       Asym_SelectFreq=glucoCESTSpectumValue_Negative-glucoCESTSpectumValue_Positive;
            %          Asym_Division=glucoCESTSpectumValue_Negative./glucoCESTSpectumValue_Positive;
                end
                SelectFreqAsymResult(k,l)=Asym_SelectFreq;
    end      
end
%  end
 toc
 cd (path);
% Base_SplineIntegralResults(Base_SplineIntegralResults>100)=0;
figure,imagesc(rot90(SelectFreqAsymResult,3),[-0.1 0.1]);colormap hot;
colorbar;
title(sprintf('SelecFreqAsymMap:%s',CESTdata.header.patientID));
% saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string')));
save('SelectFreqAsymResult.mat','SelectFreqAsymResult');
% save('BaseMap+-5ppm:spline.mat','Base_AcqData','Base_SplineCorrectedData','Base_SplineIntegralResults');
% cd ('/Users/clairezhao/brukermatlabscripts/CEST_DataAnalysis/Matlab_Analysis/InfusionGlucoCEST/');
% ImageViewerStudio(SplineIntegralResults);
end
