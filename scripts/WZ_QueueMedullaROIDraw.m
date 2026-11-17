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
             for k=1:1000
               if isfile('CESTROIValue_5ppmCtx.mat')==1
                 if isfile('spectrum.mat')==1
                        counter=counter+1;
                        attachName='_ctx';
                        dinfo = dir('spectrum.mat');
                        thisfile = dinfo.name;
                        [pathstr, base, ext] = fileparts(thisfile);
                        newfid = fullfile(pathstr, [base attachName ext]);
                        movefile(thisfile, newfid);
                        ValidCESTspectrumValue{counter,1}=SubjectPath(i).name;
                        ValidCESTspectrumValue{counter,2}=S(j).name;
                            if isfile('ROI_medulla.mat')==1
                                disp('already analyzed!')
                            elseif isfile('ROI_medulla.mat')==0
                            cd('/Users/clairezhao/Documents/MATLAB/CEST/');
                            ROI=QueueROIDraw(path);
                            CESTROI.ROI=ROI;
                            cd(path);
                            save(sprintf('ROI%s.mat',CESTROI.name),'-struct','CESTROI');
                            end
                            
                 elseif isfile('spectrum5ppmCtx.mat')==1
                        counter=counter+1;
                        replaceName='spectrum_ctx';
                        dinfo = dir('spectrum5ppmCtx.mat');
                        thisfile = dinfo.name;
                        [pathstr, base, ext] = fileparts(thisfile);
                        newfid = fullfile(pathstr, [replaceName ext]);
                        movefile(thisfile, newfid);
                        ValidCESTspectrumValue{counter,1}=SubjectPath(i).name;
                        ValidCESTspectrumValue{counter,2}=S(j).name;
                            if isfile('ROI_medulla.mat')==1
                                disp('already analyzed!')
                            elseif isfile('ROI_medulla.mat')==0
                            cd('/Users/clairezhao/Documents/MATLAB/CEST/');
                            ROI=QueueROIDraw(path);
                            CESTROI.ROI=ROI;
                            cd(path);
                            save(sprintf('ROI%s.mat',CESTROI.name),'-struct','CESTROI');
                            end
   
                 elseif isfile('spectrum.mat')==0
                     
                        counter=counter+1;
                        MissingCESTspectrumValue{counter,1}=SubjectPath(i).name;
                        MissingCESTspectrumValue{counter,2}=S(j).name;       
                            if isfile('ROI_medulla.mat')==1
                                disp('already analyzed!')
                            elseif isfile('ROI_medulla.mat')==0
                            cd('/Users/clairezhao/Documents/MATLAB/CEST/');
                            ROI=QueueROIDraw(path);
                            CESTROI.ROI=ROI;
                            cd(path);
                            save(sprintf('ROI%s.mat',CESTROI.name),'-struct','CESTROI');
                            end
                 end
               end        
              end
           end
      end
end

cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/')
save(sprintf('Results_%s.mat',datetime),'MissingCESTspectrumValue','ValidCESTspectrumValue','CESTspectrumHeader');