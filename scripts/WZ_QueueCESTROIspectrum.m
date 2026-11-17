%Wan-Ting Zhao 2021-04-29
%Document on how to pipeline my dataset
%0. Download data from cloud and store in always the same folder.
%   Run 'RenameBruker.mat', 
%   manual re-named on each measurement folder,ex: Feces, Fat, 2wk
%1. Run 'WZ_QueuePerfusioT1T2Map.mat'
%2. Run 'WZ_QueueCESTROIspectrum.mat'
%3. Run ''

%
clear all;close all;clc;
addpath('/Users/clairezhao/Documents/MATLAB/CEST/');
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/NewData/new';
cd(Basepath);
SubjectPath=dir('*MIBIO*')
% SubjectPath=dir('*KHTCH*')
for k=1:2
   if k==1 CESTROI.name=('ctx'); 
   end 
   if k==2 CESTROI.name=('medulla'); 
   end
    
    for i=1:numel(SubjectPath)
        SubjectMeasurementPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);% measurement path    
        cd(SubjectMeasurementPath);
         S = dir('*CEST_base*');%showing how many element are there in the folder     
         for j=1:numel(S)
                if S(j).isdir==1
                 path=fullfile(S(j).folder,filesep,S(j).name);
                 cd(path);
                        if isfile(sprintf('CESTROIValue_%s.mat',CESTROI.name))==1
                        load(sprintf('CESTROIValue_%s.mat',CESTROI.name))
                        CESTROI.ROI=ROI;
                        save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');

                        elseif isfile(sprintf('ROI_%s.mat',CESTROI.name))==0
                            ROI=QueueROIDraw(path);
                            CESTROI.ROI=ROI;
                            % cd(path);
                            save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
                        end
                            % Queue Spectrum_xxx.mat
                            % load(sprintf('ROI_%s.mat',CESTROI.name));
                            if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                                if isfile(sprintf('ROI_%s.mat',CESTROI.name))==1
                                    load(sprintf('ROI_%s.mat',CESTROI.name));
                                    if isfile('BaseMap+-5ppm:spline.mat')==0
                                        CESTAsymMap(path,0.5,2);
                                        WZ_ROISpectrum(path,ROI,CESTROI.name);
                                    else
                                        try load('BaseMap+-5ppm:spline.mat', 'Base_SplineCorrectedData')
                                            WZ_ROISpectrum(path,ROI,CESTROI.name);
                                        catch
                                            CESTAsymMap(path,0.5,2);
                                            WZ_ROISpectrum(path,ROI,CESTROI.name);
                                        end
                                    end

                                end
                            %     save(sprintf('spectrum_%s.mat',CESTROI.name),'-struct','CESTROI');
                            end

                end
             cd(path);
             cd ..
             S2=dir('*CEST_infu*');
             pathInfu=fullfile([S2.folder '/' S2.name]);
             cd(pathInfu);
                        if isfile(sprintf('CESTROIValue_%s.mat',CESTROI.name))==1
                        load(sprintf('CESTROIValue_%s.mat',CESTROI.name))
                        CESTROI.ROI=ROI;
                        save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');

                        elseif isfile(sprintf('ROI_%s.mat',CESTROI.name))==0
                            ROI=QueueROIDraw(pathInfu);
                            CESTROI.ROI=ROI;
                            % cd(path);
                            save(sprintf('ROI_%s.mat',CESTROI.name),'-struct','CESTROI');
                        end
                            % Queue Spectrum_xxx.mat
                            % load(sprintf('ROI_%s.mat',CESTROI.name));
                            if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==0
                                if isfile(sprintf('ROI_%s.mat',CESTROI.name))==1
                                    load(sprintf('ROI_%s.mat',CESTROI.name));
                                    if isfile('BaseMap+-5ppm:spline.mat')==0
                                        CESTAsymMap(pathInfu,0.5,2);
                                        WZ_ROISpectrum(pathInfu,ROI,CESTROI.name);
                                    else
                                        try load('BaseMap+-5ppm:spline.mat', 'Base_SplineCorrectedData')
                                            WZ_ROISpectrum(pathInfu,ROI,CESTROI.name);
                                        catch
                                            CESTAsymMap(pathInfu,0.5,2);
                                            WZ_ROISpectrum(pathInfu,ROI,CESTROI.name);
                                        end
                                    end
                                end
                            %     save(sprintf('spectrum_%s.mat',CESTROI.name),'-struct','CESTROI');
                            end

          end
    end

end