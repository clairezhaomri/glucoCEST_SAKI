%% 0.5-2ppm Queue to generate CESTMap and Store Asymmetric Spectrum integral
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*MIBIO*')
% SubjectPath=dir('*MIBIO-00725*')
StartContrastPpm=0.5;
EndContrastPpm=2;


% counter=0;
%      IntegralResultValue=cell(350,4);
%      IntegralResultHeader={'AnimalID','protocolName','ACQ_time'};
for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
              for k=1:1000
                  if isfile('BaseMap+-5ppm:spline.mat')==0
%                     counter=counter+1;  
                    cd('/Users/clairezhao/Documents/MATLAB/CEST/')
%                     try%try to go into and load the data#
                         CESTAsymMap(path,StartContrastPpm,EndContrastPpm);

%                           [AnimalID,ACQ_scan_name,ACQ_time,Base_SplineIntegralResults] = CESTAsymMap(path,StartContrastPpm,EndContrastPpm);
%                           IntegralResultValue{counter,1}=AnimalID;
%                           IntegralResultValue{counter,2}=ACQ_scan_name;
%                           IntegralResultValue{counter,3}=ACQ_time;
%                           IntegralResultValue{counter,4}=RAW_x;
%                           IntegralResultValue{counter,5}=RAW_Y;
%                           IntegralResultValue{counter,6}=xaxis4k;
%                           IntegralResultValue{counter,7}=ROIMean4K;

%                           IntegralResultValue{counter,4}=Base_SplineIntegralResults;
%                     catch% if there is missing 2dseq... 
%                         disp('No available data!')
% %                     end
%                   else
%                     disp('No available data!')
                 end
              end
           end
     end
%      close all;
end
% 
% cd('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison')
% ResultDataPath=fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm)));
% if isfolder(ResultDataPath)==0
% mkdir(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm));
% cd (fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
% else
% cd (fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
% end
% 
% save(sprintf('CESTAsymResults_%s.mat',datetime),'IntegralResultValue','IntegralResultHeader');

