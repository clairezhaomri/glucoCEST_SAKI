%% 0-2ppm input ROI_ctx and Queue spectrum_xxx.mat as output
clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData';
cd(Basepath);
SubjectPath=dir('*MIBIO*')
CESTROI.name=('ctx');


for i=1:numel(SubjectPath)
    fprintf('Iteration_counts:%3.3i/%3.3i\n', i, numel(SubjectPath)); 
    SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
    cd(SubjectMeasurementPath);
     S = dir('*CEST*');%showing how many element are there in the folder     
     for j=1:numel(S) 
           if S(j).isdir==1
             path=fullfile(S(j).folder,filesep,S(j).name);
             cd(path);
%              if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==1
%                 disp('already analyzed');
%                 continue;
             if isfile(sprintf('ROI_%s.mat',CESTROI.name))==0
                 continue;
             else
                        if isfile('BaseMap+-5ppm:spline.mat')==1
                            load('BaseMap+-5ppm:spline.mat');
                             if (exist ('Base_AcqData','var'))==1
                                 cd(path)
                                 load(sprintf('ROI_%s.mat',CESTROI.name));
                                 cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                  WZ_ROISpectrum(path,ROI,CESTROI.name);
                             elseif (exist ('Base_AcqData','var'))==0
                                 cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                 Data=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']); 
                                    if size(Data.dims(1))==size(Data.dims(2))
                                        try
                                            CESTAsymMap(path,0.5,2);
                                        catch
                                            disp('dimension mismatch!');
                                        end
                                    end
                                 load(sprintf('ROI_%s.mat',CESTROI.name));
                                 cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                         WZ_ROISpectrum(path,ROI,CESTROI.name);
                             end
                        else
                        end


                        if isfile('CESTBaseMap(spline).mat')==1
                            if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==1
                                continue;
                            elseif isfile('BaseMap+-5ppm:spline.mat')==1
                                    delete('CESTBaseMap(spline).mat');
                                    continue;
                            else
                                    load ('CESTBaseMap(spline).mat');

                                     if (exist ('Base_AcqData','var'))==1
                                         cd(path);
                                         load(sprintf('ROI_%s.mat',CESTROI.name));
                                         cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                         WZ_ROISpectrum(path,ROI,CESTROI.name);
                                     elseif (exist ('Base_AcqData','var'))==0
                                         cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                         Data=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']); 
                                            if size(Data.dims(1))==size(Data.dims(2))
                                                try
                                                    CESTAsymMap(path,0.5,2);
                                                catch
                                                    disp('dimension mismatch!');
                                                end
                                            end
                                         load(sprintf('ROI_%s.mat',CESTROI.name));
                                         cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                         WZ_ROISpectrum(path,ROI,CESTROI.name);
                                         delete('CESTBaseMap(spline).mat');
                                     end
                            end
                            
                        end
                          
                        if isfile('CESTBaseMap(pchip).mat')==1
                             if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==1
                                continue;
                            elseif isfile('BaseMap+-5ppm:spline.mat')==1
                                    delete('CESTBaseMap(pchip).mat');
                                    continue;
                             else
                                 cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                 Data=MRIDataLoader.LoadFrom2DSeq([path filesep 'pdata' filesep '1' filesep '2dseq']); 
                                    if size(Data.dims(1))==size(Data.dims(2))
                                        try
                                            CESTAsymMap(path,0.5,2);
                                        catch
                                            disp('dimension mismatch!');
                                        end
                                    end
                                 cd(path)
                                 load(sprintf('ROI_%s.mat',CESTROI.name));
                                 cd('/Users/clairezhao/Documents/MATLAB/CEST');
                                         WZ_ROISpectrum(path,ROI,CESTROI.name);
                                 delete('CESTBaseMap(pchip).mat');
                             end
                        end
             end
           end
     end
end