clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('KHTCH-11043');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=0.4;



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
                 if isfile('spectrum.mat')==1
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/Users/clairezhao/Documents/MATLAB/CEST')
%                     try%try to go into and load the data#
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm);
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

cd('/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');


%%
clear all;close all;clc;
Basepath='/data/projects/bruker-hoerr-cest/data/NewData/';
cd(Basepath);
SubjectPath=dir('*MIBIO*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=0.66;



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
                 if isfile('spectrum.mat')==1
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/data/projects/bruker-hoerr-cest/scripts/wzhao/CEST_DataAnalysis/DataPipelineQueue')
%                     try%try to go into and load the data#
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm);
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

cd('/data/projects/bruker-hoerr-cest/documents/')
ResultDataPath=fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');

%%

clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*MIBIO*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=1.28;



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
                 if isfile('spectrum5ppmCtx.mat')==1
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/Users/clairezhao/Documents/MATLAB/CEST')
%                     try%try to go into and load the data#
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm);
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

cd('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison/CtrlInfection/test/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison/CtrlInfection/test/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/acuteComparison/CtrlInfection/test/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');

%%
clear all;close all;clc;
Basepath='/data/projects/bruker-hoerr-cest/data/NewData/';
cd(Basepath);
SubjectPath=dir('*NutritionCEST*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=2.08;



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
                 if isfile('spectrum.mat')==1
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/data/projects/bruker-hoerr-cest/scripts/wzhao/CEST_DataAnalysis/DataPipelineQueue')
%                     try%try to go into and load the data#
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm);
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

cd('/data/projects/bruker-hoerr-cest/documents/')
ResultDataPath=fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');

%%
clear all;close all;clc;
Basepath='/data/projects/bruker-hoerr-cest/data/NewData/';
cd(Basepath);
SubjectPath=dir('*NutritionCEST*');
% SubjectPath=dir('*MIBIO-00730*')
glucoCESTppm=2.88;



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
                 if isfile('spectrum.mat')==1
                   counter=counter+1;  
                    CESTspectrumResultValue{counter,1}=SubjectPath(i).name;
                    CESTspectrumResultValue{counter,2}=S(j).name;
                    cd('/data/projects/bruker-hoerr-cest/scripts/wzhao/CEST_DataAnalysis/DataPipelineQueue')
%                     try%try to go into and load the data#
                        [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz(path,glucoCESTppm);
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

cd('/data/projects/bruker-hoerr-cest/documents/')
ResultDataPath=fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('AsymGlucoCESTppm_%3.2f',glucoCESTppm));
else
cd (fullfile('/data/projects/bruker-hoerr-cest/documents/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
end
save(sprintf('CESTAsymResults_%3.2fppm%s.mat',glucoCESTppm,datetime),'CESTspectrumResultValue','spectrumResultHeader');
