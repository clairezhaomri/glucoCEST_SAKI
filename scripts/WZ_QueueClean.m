%% 0.5-2ppm Queue to generate CESTMap and Store Asymmetric Spectrum integral
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
% SubjectPath=dir('*SEK00*');



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
                  if isfile('BaseMap+-5ppm:spline.mat')==1
                    delete('BaseMap+-5ppm:spline.mat')
                  end
                    if isfile('BaseMap+-5ppm_spline.mat')==1
                    delete('BaseMap+-5ppm_spline.mat')
                    end

                 
              
           end
     end
%      close all;
end
disp('done');
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

