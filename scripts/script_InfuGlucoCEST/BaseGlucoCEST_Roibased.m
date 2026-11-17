%% load data

% CEST part
CESTpath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/20200319_211606_MIBIO_00534_1_1FatAcute/7CESTbase/';
CESTdatapath=fullfile(CESTpath,filesep,'pdata',filesep,'1');
CESTdata=MRIDataLoader.LoadFrom2DSeq(CESTdatapath); 
CEST_Offsets=CESTdata.headerHistory{1}.method.MT_Offsets_NoM0;
CESTdata.Permute('read','line','echo');
handles.CEST.CEST=CESTdata;
handles.CEST.DataSlope = CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope;


% %if no S0data
%  S0=CESTdata.SubVolume('echo',);
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


%% manual draw Roi-based analysis
% imagesc(squeeze(CESTdata(:,:,1)));
% ROI=roipoly;

[xSort, xIndex] = sort(CEST_Offsets);
locscaleddata=handles.CEST.CEST.vol ./ localCEST_S0_Vol;
locscaleddata=locscaleddata(:,:,xIndex);
    
RoiMean=zeros(handles.CEST.CEST.headerHistory{1}.method.Number_Offsets,1);
RoiStd=zeros(handles.CEST.CEST.headerHistory{1}.method.Number_Offsets,1);



for echoi = 1:handles.CEST.CEST.headerHistory{1}.method.Number_Offsets
   singleecho=locscaleddata(:,:,echoi);
   RoiMean(echoi)=mean(singleecho(ROI_kidney));
   RoiStd(echoi)=std(singleecho(ROI_kidney));
end


figure,errorbar(xSort,RoiMean,RoiStd,'-*');

pause(1);