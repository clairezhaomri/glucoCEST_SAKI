%% MultipleQueue
clear all;close all;clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*20*');
CESTROI.name=('ctx');


counter=0;
     IntegralResultValue=cell(200,9);
     IntegralResultHeader={'AnimalID','phase','diet','CESTBaseAPT3.5ppm','CESTAPTInfu3.5ppm','glucoCEST3.5ppm','CESTNOEBase-3.5ppm','CESTNOEInfu-3.5ppm','glucoCEST-3.5ppm'};
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

                          name=CESTROI.name;
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(path,3.5,name,'APT');                         
                          IntegralResultValue{counter,1}=AnimalID;
                          IntegralResultValue{counter,4}=Asym_SelectFreq;
                          
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(path,-3.5,name,'NOE');                          
                          IntegralResultValue{counter,7}=Asym_SelectFreq;
                      
                      cd(path);
                      cd ..
                      S2=dir('*CEST_infu*');
                      if numel(S2)==0  
                      IntegralResultValue{counter,5}='x';
                      IntegralResultValue{counter,6}='x';
                      IntegralResultValue{counter,8}='x';
                      IntegralResultValue{counter,9}='x';
                          continue;
                      else
                      pathInfu=fullfile([S2.folder '/' S2.name]);
                      cd(pathInfu);
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(pathInfu,3.5,name,'APT');                          
                              IntegralResultValue{counter,5}=Asym_SelectFreq;
                              IntegralResultValue{counter,6}=IntegralResultValue{counter,5}-IntegralResultValue{counter,4};
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(pathInfu,-3.5,name,'NOE');                          
                              IntegralResultValue{counter,8}=Asym_SelectFreq;
                              IntegralResultValue{counter,9}=IntegralResultValue{counter,8}-IntegralResultValue{counter,7};
                      end
                    
                  end
              end
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

%% MultipleQueue
clear all;close all;clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*20*');
CESTROI.name=('medulla');


counter=0;
     IntegralResultValue=cell(200,9);
     IntegralResultHeader={'AnimalID','phase','diet','CESTBaseAPT3.5ppm','CESTAPTInfu3.5ppm','glucoCEST3.5ppm','CESTNOEBase-3.5ppm','CESTNOEInfu-3.5ppm','glucoCEST-3.5ppm'};
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

                          name=CESTROI.name;
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(path,3.5,name,'APT');                         
                          IntegralResultValue{counter,1}=AnimalID;
                          IntegralResultValue{counter,4}=Asym_SelectFreq;
                          
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(path,-3.5,name,'NOE');                          
                          IntegralResultValue{counter,7}=Asym_SelectFreq;
                      
                      cd(path);
                      cd ..
                      S2=dir('*CEST_infu*');
                      if numel(S2)==0  
                      IntegralResultValue{counter,5}='x';
                      IntegralResultValue{counter,6}='x';
                      IntegralResultValue{counter,8}='x';
                      IntegralResultValue{counter,9}='x';
                          continue;
                      else
                      pathInfu=fullfile([S2.folder '/' S2.name]);
                      cd(pathInfu);
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(pathInfu,3.5,name,'APT');                          
                              IntegralResultValue{counter,5}=Asym_SelectFreq;
                              IntegralResultValue{counter,6}=IntegralResultValue{counter,5}-IntegralResultValue{counter,4};
                          [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(pathInfu,-3.5,name,'NOE');                          
                              IntegralResultValue{counter,8}=Asym_SelectFreq;
                              IntegralResultValue{counter,9}=IntegralResultValue{counter,8}-IntegralResultValue{counter,7};
                      end
                    
                  end
              end
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