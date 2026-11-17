%% MultipleQueue
clear all;close all;clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Volumes/ClaireData/SophieWobst/';
cd(Basepath);
% SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
SubjectPath=[dir('*SEK*')];

CESTROI.name=('ctx');

counter=0;
     IntegralResultValue=cell(30,13);
     IntegralResultHeader={'AnimalID','experimentType','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
         'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
         'rNOEBase-3.5ppm'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*bCEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%               for k=1:1000
%                   if isfile(sprintf('CESTROIValue_%s.mat',CESTROI.name))==0
%                         continue;
%                   else
                    counter=counter+1;
                        %write phase
                        header=[S(j).folder];
                        if contains(header,'CLP')
                           IntegralResultValue{counter,2}='CLP';
                        else 
                           IntegralResultValue{counter,2}='Ctrl';
                        end
                       
%                     cd('/Users/clairezhao/Documents/MATLAB/CEST');
%                     try%try to go into and load the data#
                          name=CESTROI.name;
                         [AnimalID,IntegralValue] = SW_CESTAsymSplineIntegral(path,0.5,2,name);
                          IntegralResultValue{counter,1}=AnimalID;
                          IntegralResultValue{counter,3}=IntegralValue;
                          
                          [AnimalID,IntegralValue] = SW_CESTAsymSplineIntegral(path,1.2,2.2,name);
                          IntegralResultValue{counter,6}=IntegralValue;
                          
                          [Asym_SelectFreq]=SW_CESTSpectrumInSelectedHz(path,1.28,name);
                          IntegralResultValue{counter,9}=Asym_SelectFreq;
                          [Asym_SelectFreq]=SW_CESTSpectrumInSelectedHz(path,3.5,name);
                          IntegralResultValue{counter,12}=Asym_SelectFreq;
                          [Asym_SelectFreq]=SW_CESTSpectrumInSelectedHz(path,-3.5,name);
                          IntegralResultValue{counter,13}=Asym_SelectFreq;
                      
                          
%                       cd(path);
%                       cd ..
%                       S2=dir('*iCEST*');
%                       if numel(S2)==0  
%                       IntegralResultValue{counter,4}='x';
%                       IntegralResultValue{counter,5}='x';
%                       IntegralResultValue{counter,7}='x';
%                       IntegralResultValue{counter,8}='x';
%                       IntegralResultValue{counter,10}='x';
%                       IntegralResultValue{counter,11}='x';
% 
%                       
%                           continue;
%                       else
%                       pathInfu=fullfile([S2.folder '/' S2.name]);
%                       cd(pathInfu);
%                       [AnimalID,IntegralValue] = SW_CESTAsymSplineIntegral(pathInfu,0.5,2,name);
%                       IntegralResultValue{counter,4}=IntegralValue;
%                       IntegralResultValue{counter,5}=IntegralResultValue{counter,4}-IntegralResultValue{counter,3};
%                       [AnimalID,IntegralValue] = SW_CESTAsymSplineIntegral(pathInfu,1.2,2.2,name);
%                       IntegralResultValue{counter,7}=IntegralValue;
%                       IntegralResultValue{counter,8}=IntegralResultValue{counter,7}-IntegralResultValue{counter,6};
%                       
%                       [Asym_SelectFreq]=SW_CESTSpectrumInSelectedHz(pathInfu,1.28,name);
%                       IntegralResultValue{counter,10}=Asym_SelectFreq;
%                       IntegralResultValue{counter,11}=IntegralResultValue{counter,10}-IntegralResultValue{counter,9};

%                       end
                    
%                   end
%               end
           end
     end
end
%      close all;


cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTResult%s',date)));
if isfolder(ResultDataPath)==0
mkdir(sprintf('glucoCESTResult%s',date));
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTResult%s',date))));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTResult%s',date))));
end

save('SW-CESTAsymResults.mat','IntegralResultValue','IntegralResultHeader');

