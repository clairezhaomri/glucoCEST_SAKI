%% MultipleQueue
clear all;close all;clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/';
cd(Basepath);
SubjectPath=dir('*20200721_233343_Phantom_glucose_multi_1_1*');
CESTROI.name=('*phantom*');


counter=0;
     IntegralResultValue=cell(200,23);
     IntegralResultHeader={'IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
         'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_base*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
              for k=1:1000
                  if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                        continue;
                  else
                    counter=counter+1;

                          name=CESTROI.name;
                         [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(path,0.5,2,name);
                          IntegralResultValue{counter,1}=AnimalID;
%                           IntegralResultValue{counter,5}=RAW_x;
%                           IntegralResultValue{counter,6}=RAW_Y;
                          IntegralResultValue{counter,9}=IntegralValue;
%                           
                          [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(path,1.2,2.2,name);
                          IntegralResultValue{counter,12}=IntegralValue;
                          [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq]=WZ_CESTSpectrumInSelectedHz(path,1.28,name);
                          IntegralResultValue{counter,15}=Asym_SelectFreq;
                          [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq]=WZ_CESTSpectrumInSelectedHz(path,3.5,name);
                          IntegralResultValue{counter,18}=Asym_SelectFreq;
                          [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq]=WZ_CESTSpectrumInSelectedHz(path,-3.5,name);
                          IntegralResultValue{counter,21}=Asym_SelectFreq;
                      
%                       cd(path);
%                       cd ..
%                       S2=dir('*CEST_infu*');
%                       if numel(S2)==0  
%                       IntegralResultValue{counter,7}='x';
%                       IntegralResultValue{counter,8}='x';
%                       IntegralResultValue{counter,10}='x';
%                       IntegralResultValue{counter,11}='x';
%                       IntegralResultValue{counter,13}='x';
%                       IntegralResultValue{counter,14}='x';
%                       IntegralResultValue{counter,16}='x';
%                       IntegralResultValue{counter,17}='x';
%                       
%                           continue;
%                       else
%                       pathInfu=fullfile([S2.folder '/' S2.name]);
%                       cd(pathInfu);
%                       [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(pathInfu,0.5,2,name);
%                       IntegralResultValue{counter,7}=RAW_x;
%                       IntegralResultValue{counter,8}=RAW_Y;
%                       IntegralResultValue{counter,10}=IntegralValue;
%                       IntegralResultValue{counter,11}=IntegralResultValue{counter,10}-IntegralResultValue{counter,9};
%                       [AnimalID,ACQ_scan_name,ACQ_time,RAW_x,RAW_Y,xaxis4k,ROIMean4K,IntegralValue] = WZ_CESTAsymSplineIntegral(pathInfu,1.2,2.4,name);
%                       IntegralResultValue{counter,13}=IntegralValue;
%                       IntegralResultValue{counter,14}=IntegralResultValue{counter,13}-IntegralResultValue{counter,12};
%                       [glucoCESTSpectumValue_Positive,glucoCESTSpectumValue_Negative,Asym_SelectFreq]=WZ_CESTSpectrumInSelectedHz(pathInfu,1.28,name);
%                       IntegralResultValue{counter,16}=Asym_SelectFreq;
%                       IntegralResultValue{counter,17}=IntegralResultValue{counter,16}-IntegralResultValue{counter,15};

%                       end
                    
                  end
              end
           end
     end
end
%      close all;


cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('phantomCESTResult%s',date)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('glucoCESTResult%s',date));
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('phantomCESTResult%s',date))));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('phantomCESTResult%s',date))));
end

save(sprintf('CESTAsymResults_%s%s.mat',CESTROI.name,datetime),'IntegralResultValue','IntegralResultHeader');