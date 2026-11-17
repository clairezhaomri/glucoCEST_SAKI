


%% 1.2-2.4ppm Original script
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new';
cd(Basepath);
SubjectPath=dir('*MIBIO*')
% SubjectPath2=dir('*KHTCH*')

% SubjectPath=dir('*MIBIO-00725*')
StartContrastPpm=1.2;
EndContrastPpm=2.4;
CESTROI.name=('medulla');


counter=0;
     IntegralResultValue=cell(150,8);
     IntegralResultHeader={'AnimalID','protocolName','ACQ_time','AcquiredOffsets','AcquiredMeanValue','SplinedOffsets_x','SplinedMeanValue_y','IntegralValue'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
              for k=1:1000
                  if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                        continue;
                  else
                    counter=counter+1;  

                    cd('/Users/clairezhao/Documents/MATLAB/CEST');
%                     try%try to go into and load the data#
                          name=CESTROI.name;
                         [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(path,StartContrastPpm,EndContrastPpm,name);
                          IntegralResultValue{counter,1}=AnimalID;
                          IntegralResultValue{counter,2}=ACQ_scan_name;
                          IntegralResultValue{counter,3}=ACQ_time;
                          IntegralResultValue{counter,4}=RAW_x;
                          IntegralResultValue{counter,5}=RAW_Y;
                          IntegralResultValue{counter,6}=xaxis4k;
                          IntegralResultValue{counter,7}=ROIMean4K;
                          IntegralResultValue{counter,8}=IntegralValue;
%                     catch% if there is missing 2dseq... 
%                         disp('No available data!')
% %                     end
%                   else
%                     disp('No available data!')
                  end
              end
           end
     end
end
%      close all;


cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm));
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
end

save(sprintf('CESTAsymResults_%s%s.mat',CESTROI.name,datetime),'IntegralResultValue','IntegralResultHeader');


%% 0.5-2ppm Original script
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new';
cd(Basepath);
SubjectPath=dir('*KHTCH*')
% SubjectPath2=dir('*KHTCH*')

% SubjectPath=dir('*MIBIO-00725*')
StartContrastPpm=0.5;
EndContrastPpm=2;
CESTROI.name=('ctx');


counter=0;
     IntegralResultValue=cell(150,8);
     IntegralResultHeader={'AnimalID','protocolName','ACQ_time','AcquiredOffsets','AcquiredMeanValue','SplinedOffsets_x','SplinedMeanValue_y','IntegralValue'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
              for k=1:1000
                  if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                        continue;
                  else
                    counter=counter+1;  

                    cd('/Users/clairezhao/Documents/MATLAB/CEST');
%                     try%try to go into and load the data#
                          name=CESTROI.name;
                         [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(path,StartContrastPpm,EndContrastPpm,name);
                          IntegralResultValue{counter,1}=AnimalID;
                          IntegralResultValue{counter,2}=ACQ_scan_name;
                          IntegralResultValue{counter,3}=ACQ_time;
                          IntegralResultValue{counter,4}=RAW_x;
                          IntegralResultValue{counter,5}=RAW_Y;
                          IntegralResultValue{counter,6}=xaxis4k;
                          IntegralResultValue{counter,7}=ROIMean4K;
                          IntegralResultValue{counter,8}=IntegralValue;
%                     catch% if there is missing 2dseq... 
%                         disp('No available data!')
% %                     end
%                   else
%                     disp('No available data!')
                  end
              end
           end
     end
end
%      close all;


cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm));
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTIntegral@%3.1f-%3.1fppm',StartContrastPpm,EndContrastPpm))));
end

save(sprintf('CESTAsymResults_%s%s.mat',CESTROI.name,datetime),'IntegralResultValue','IntegralResultHeader');



%%

clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*MIBIO*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=1.28;
CESTROI.name=('medulla');



     counter=0;
     CESTspectrumResultValue=cell(200,6);
     spectrumResultHeader={'AnimalID','protocolName','glucoCESTppm','glucoCESTSpectumValue_Positive','glucoCESTSpectumValue_Negative','Asym_SelectFreq'};
for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
             for k=1:1000
                 if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                     continue;
                 else
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/Users/clairezhao/Documents/MATLAB/CEST')
%                     try%try to go into and load the data#
                         name=CESTROI.name;
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm,name);
                        CESTspectrumResultValue{counter,3}=glucoCESTppm;
                        CESTspectrumResultValue{counter,4}=glucoCESTSpectumValue_Positive;
                        CESTspectrumResultValue{counter,5}=glucoCESTSpectumValue_Negative;
                        CESTspectrumResultValue{counter,6}=Asym_SelectFreq;

%                     catch% if there is missing 2dseq... 
%                         disp('No available data!')
%                     end        
                 end
             end
           else
                    disp('No available data!')
           
           end
     end
%      close all;
end

cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');


%%
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*KHTCH*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=1.28;
CESTROI.name=('medulla');



     counter=0;
     CESTspectrumResultValue=cell(200,6);
     spectrumResultHeader={'AnimalID','protocolName','glucoCESTppm','glucoCESTSpectumValue_Positive','glucoCESTSpectumValue_Negative','Asym_SelectFreq'};
for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
             for k=1:1000
                 if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                     continue;
                 else
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/Users/clairezhao/Documents/MATLAB/CEST')
%                     try%try to go into and load the data#
                         name=CESTROI.name;
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm,name);
                        CESTspectrumResultValue{counter,3}=glucoCESTppm;
                        CESTspectrumResultValue{counter,4}=glucoCESTSpectumValue_Positive;
                        CESTspectrumResultValue{counter,5}=glucoCESTSpectumValue_Negative;
                        CESTspectrumResultValue{counter,6}=Asym_SelectFreq;

%                     catch% if there is missing 2dseq... 
%                         disp('No available data!')
%                     end        
                 end
             end
           else
                    disp('No available data!')
           
           end
     end
%      close all;
end

cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');

