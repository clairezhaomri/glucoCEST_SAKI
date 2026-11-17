%% MultipleQueue
clear all;close all;
% clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
%  SubjectPath=dir('*13C*');

CESTROI.name=('ctx');
%%
counter=0;
     IntegralResultValue=cell(200,23);
     IntegralResultHeader={'AnimalID','phase','diet','experimentType','AcquiredBaseOffsets','AcquiredBaseROIMeanValue','AcquiredInfuOffsets','AcquiredInfuROIMeanValue','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
         'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
         'rNOEBase-3.5ppm'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_base*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%               for k=1:1000
%                   if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
%                         continue;
%                   else
                    counter=counter+1;
                        %write phase
                        header=[S(j).folder];
                        if contains(header,'Acute')
                           IntegralResultValue{counter,2}='acute';
                        elseif contains(header,'2wk')
                           IntegralResultValue{counter,2}='2wk';
                        elseif contains(header,'8wk')
                           IntegralResultValue{counter,2}='8wk';
                        end
                        
                        header=[S(j).folder];
                        if contains(header,'Fat')
                           IntegralResultValue{counter,3}='Fat';
                        elseif contains(header,'Protein')
                           IntegralResultValue{counter,3}='Protein';
                        elseif contains(header,'Normal')
                           IntegralResultValue{counter,3}='Normal';
                        end

                        
                        if contains(header,'Feces')
                           IntegralResultValue{counter,4}='Feces Infect';
                        elseif contains(header,'Ctrl')
                           IntegralResultValue{counter,4}='Ctrl Infect';
                        end
%                     cd('/Users/clairezhao/Documents/MATLAB/CEST');
%                     try%try to go into and load the data#
                           name=CESTROI.name;
%                           [AnimalID,IntegralValue,RAW_x,RAW_Y] = WT_CESTAsymSplineIntegral(path,0,1,name);%ESTbase@0.5-2ppm

%                           [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,0,1,name);%ESTbase@0.5-2ppm
                          [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,2,3.6,name);%ESTbase@0.5-2ppm

                          IntegralResultValue{counter,1}=AnimalID;
%                           IntegralResultValue{counter,5}=RAW_x;
%                           IntegralResultValue{counter,6}=RAW_Y; %if an
%                           %outout is needed, corresponding change in
%                           function are reqired.
                          IntegralResultValue{counter,9}=IntegralValue;
                          
%                           [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,0.2,0.8,name);%lactate
%                           IntegralResultValue{counter,12}=IntegralValue;
                          
                          [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,3,name);
                          IntegralResultValue{counter,15}=Asym_SelectFreq;
%                           [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,0.5,name);%;lactate
%                           IntegralResultValue{counter,18}=Asym_SelectFreq;
%                           [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,0.5,name);
%                           IntegralResultValue{counter,21}=Asym_SelectFreq;
                      
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
%                       [AnimalID,IntegralValue,RAW_x,RAW_Y] = WT_CESTAsymSplineIntegral(pathInfu,0,1,name);
% %                       IntegralResultValue{counter,7}=RAW_x;
% %                       IntegralResultValue{counter,8}=RAW_Y;
%                       IntegralResultValue{counter,10}=IntegralValue;
%                       IntegralResultValue{counter,11}=IntegralResultValue{counter,10}-IntegralResultValue{counter,9};
%                       [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(pathInfu,0.2,0.8,name);
%                       IntegralResultValue{counter,13}=IntegralValue;
%                       IntegralResultValue{counter,14}=IntegralResultValue{counter,13}-IntegralResultValue{counter,12};
%                       [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(pathInfu,0.4,name);
%                       IntegralResultValue{counter,16}=Asym_SelectFreq;
%                       IntegralResultValue{counter,17}=IntegralResultValue{counter,16}-IntegralResultValue{counter,15};

%                       end
                    
                  
              
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

save(sprintf('CESTAsymResultsLactate_%s%s.mat',CESTROI.name,datetime),'IntegralResultValue','IntegralResultHeader');




%% MultipleQueue
clear all;close all;
% clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=[dir('*MIBIO*'); dir('*KHTCH*')];
% SubjectPath=dir('*13C*');

CESTROI.name=('medulla');
%%
counter=0;
     IntegralResultValue=cell(200,23);
     IntegralResultHeader={'AnimalID','phase','diet','experimentType','AcquiredBaseOffsets','AcquiredBaseROIMeanValue','AcquiredInfuOffsets','AcquiredInfuROIMeanValue','IntegralCESTBase0.5-2ppm','IntegralCESTInfu0.5-2ppm',...
         'glucoCEST0.5-2ppm','IntegralCESTBase1.2-2.2ppm','IntegralCESTInfu1.2-2.2ppm','glucoCEST1.2-2.2ppm','CESTBase1.28ppm','CESTInfu1.28ppm','glucoCEST1.28ppm','APTBase3.5ppm',...
         'rNOEBase-3.5ppm'};
for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST_base*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%               for k=1:1000
%                   if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
%                         continue;
%                   else
                    counter=counter+1;
                        %write phase
                        header=[S(j).folder];
                        if contains(header,'Acute')
                           IntegralResultValue{counter,2}='acute';
                        elseif contains(header,'2wk')
                           IntegralResultValue{counter,2}='2wk';
                        elseif contains(header,'8wk')
                           IntegralResultValue{counter,2}='8wk';
                        end
                        
                        header=[S(j).folder];
                        if contains(header,'Fat')
                           IntegralResultValue{counter,3}='Fat';
                        elseif contains(header,'Protein')
                           IntegralResultValue{counter,3}='Protein';
                        elseif contains(header,'Normal')
                           IntegralResultValue{counter,3}='Normal';
                        end

                        
                        if contains(header,'Feces')
                           IntegralResultValue{counter,4}='Feces Infect';
                        elseif contains(header,'Ctrl')
                           IntegralResultValue{counter,4}='Ctrl Infect';
                        end
%                     cd('/Users/clairezhao/Documents/MATLAB/CEST');
%                     try%try to go into and load the data#
                           name=CESTROI.name;
%                           [AnimalID,IntegralValue,RAW_x,RAW_Y] = WT_CESTAsymSplineIntegral(path,0.5,2,name);%ESTbase@0.5-2ppm
%                           IntegralResultValue{counter,1}=AnimalID;
%                           IntegralResultValue{counter,5}=RAW_x;
%                           IntegralResultValue{counter,6}=RAW_Y;
%                           IntegralResultValue{counter,9}=IntegralValue;
                          
                          [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(path,2,3.6,name);
                          IntegralResultValue{counter,1}=AnimalID;

                          IntegralResultValue{counter,12}=IntegralValue;
                          
%                           [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,1.28,name);
%                           IntegralResultValue{counter,15}=Asym_SelectFreq;
                          [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,3,name);
                          IntegralResultValue{counter,18}=Asym_SelectFreq;
%                           [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(path,-3.5,name);
%                           IntegralResultValue{counter,21}=Asym_SelectFreq;
                      
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
%                       [AnimalID,IntegralValue,RAW_x,RAW_Y] = WT_CESTAsymSplineIntegral(pathInfu,2,3.6,name);
% %                       IntegralResultValue{counter,7}=RAW_x;
% %                       IntegralResultValue{counter,8}=RAW_Y;
%                       IntegralResultValue{counter,10}=IntegralValue;
% %                       IntegralResultValue{counter,11}=IntegralResultValue{counter,10}-IntegralResultValue{counter,9};
% %                       [AnimalID,IntegralValue] = WT_CESTAsymSplineIntegral(pathInfu,1.2,2.2,name);
% %                       IntegralResultValue{counter,13}=IntegralValue;
% %                       IntegralResultValue{counter,14}=IntegralResultValue{counter,13}-IntegralResultValue{counter,12};
%                       [Asym_SelectFreq]=WT_CESTSpectrumInSelectedHz(pathInfu,1.28,name);
%                       IntegralResultValue{counter,16}=Asym_SelectFreq;
% %                       IntegralResultValue{counter,17}=IntegralResultValue{counter,16}-IntegralResultValue{counter,15};

%                       end
                    
                  
              
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

save(sprintf('CESTAsymResults_%s%s.mat',CESTROI.name,datetime),'IntegralResultValue','IntegralResultHeader');

















