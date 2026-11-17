clear all, close all, clc,
%% load data
% BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/ANEST_71983__HemeCEST_07_Heme__2019-10-08_05.07.44';
% cd(BasePath);


% CEST part
CESTInfupath='/Users/clairezhao/Documents/data/pv6/CEST/SEK00-72863__CEST_Liver__2019-08-20_10.48.00/E016_wz_RC_FS_1.8uT_TrigerOutInclude_E16';
CESTInfudatapath=fullfile(CESTInfupath,filesep,'pdata',filesep,'1');
CESTInfudata=MRIDataLoader.LoadFrom2DSeq(CESTInfudatapath); 
CEST_Offsets=CESTInfudata.headerHistory{1}.method.MT_Offsets_NoM0;
CESTInfudata.Permute('read','line','echo');
handles.CEST.CESTInfu=CESTInfudata;
% handles.CEST.InfuDataSlope = CESTInfudata.headerHistory{1}.visupars.VisuCoreDataSlope;


%load highres data
T2hiresPath=dir('*T2RARE_highres_axialMS*');
T2hiresDataPath=fullfile(BasePath,T2hiresPath.name,filesep,'pdata',filesep,'1');
Hiresdata=MRIDataLoader.LoadFrom2DSeq(T2hiresDataPath); 
Hiresdata.Permute('read','line','echo');

% 
% 
% %if no S0data
% S0=CESTdata.SubVolume('echo',1);
% DataSlope_S0=CESTdata.headerHistory{1}.visupars.VisuCoreDataSlope(1);


%normalized_CESTdata=CESTdata./S0;
% % get separate S0 measurement 
% S0path='/data/projects/bruker-hoerr-cest/data/phantom_glucose_multiConcentr__CEST_multi_04__2019-02-06_17:53:37/wz_RARE_CEST_S0_E6';
% S0datapath=fullfile(S0path,filesep,'pdata',filesep,'1');
% S0data=MRIDataLoader.LoadFrom2DSeq(S0datapath);
% S0data.Permute('echo','read','line');
% handles.CEST.CEST_S0=S0data;


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
%% manual draw Roi-based analysis
I=figure,imagesc(Hiresdata(:,:,8));
tempROI=roipoly;
waitfor(I);
tempROI=tempROI(1:4:end,1:4:end);


[xSort, xIndex] = sort(CEST_Offsets);
locscaleddata=handles.CEST.CESTInfu.vol ./ localInfuCEST_S0_Vol;
locscaleddata=locscaleddata(:,:,xIndex);


RoiMeanInfu=zeros(handles.CEST.CESTInfu.headerHistory{1}.method.Number_Offsets,1);
RoiStdInfu=zeros(handles.CEST.CESTInfu.headerHistory{1}.method.Number_Offsets,1);



for echoi = 1:handles.CEST.CESTInfu.headerHistory{1}.method.Number_Offsets
   singleecho=locscaleddata(:,:,echoi);
   RoiMeanInfu(echoi)=mean(singleecho(tempROI));
   RoiStdInfu(echoi)=std(singleecho(tempROI));
end

VoxelNum=length(find(tempROI>0));
ROISemBase=ROIstdBase./VoxelNum;

 %make symmetric
       [min_value,min_idx]=min(RoiMeanInfu);
        a=xaxis(1,min_idx);
        xaxis_corrected1=xaxis-a;
        min_value=[];
        min_idx=[];
        a=[];
        [min_value,min_idx]=min(ROIMeanInfu);
        a=xaxis(1,min_idx);
        xaxis_corrected2=xaxis-a;


errorbar(xSort,RoiMeanInfu,RoiStdInfu,'r');

pause(1);