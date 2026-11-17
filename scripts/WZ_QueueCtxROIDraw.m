clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/';
cd(Basepath);
SubjectPath=dir('*MIBIO*');
% SubjectPath=dir('*MIBIO-00730*')
CESTROI.name=('ctx');
                      
     counter=0;
     MissingCESTspectrumValue=cell(200,2);
     CESTspectrumHeader={'AnimalID','protocolName'};
     ValidCESTspectrumValue=cell(200,2);
for i=1:numel(SubjectPath)
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S)         
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
                    cd(path);
%              for k=1:1000
%                 if isfile('CESTBaseMap(pchip).mat')==0
%                     continue;
%                 elseif  isfile('BaseMap+-5ppm:spline.mat')==0
%                     continue;
%                 elseif  isfile('CESTBaseMap(spline).mat')==0
%                     continue;    
%                 else    
                 if  isfile('BaseMap+-5ppm:spline.mat')==1

                            load('CESTROIValue_5ppmCtx.mat', 'ROI');
                            CESTROI.ROI=ROI;
                            cd(path);
                            save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');

                elseif isfile('CESTROIValue_5ppmCtx.mat')==0
                            if isfile('ROI_ctx.mat')==1
                                disp('already draw!')
                            elseif isfile('ROI_ctx.mat')==0
                            cd('/Users/clairezhao/Documents/MATLAB/CEST/');
                            ROI=QueueROIDraw(path);
                            CESTROI.ROI=ROI;
                            cd(path);
                            save(sprintf('ROI%s.mat',CESTROI.name),'-struct','CESTROI');
                            end
 
                else     
                     disp('no proper data')
                end
           end        
     end      
end

cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
save(sprintf('Results_%s.mat',datetime),'MissingCESTspectrumValue','ValidCESTspectrumValue','CESTspectrumHeader');