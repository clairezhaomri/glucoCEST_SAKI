%Wan-Ting Zhao Universitaetsklinikum Jena; version 01.May. 2020. 
%This script is designed for loading the existing fit results of FAIR perfusion data, and
%plot pixelwised fit result from our ROI. 


clear all;close all;clc;
Basepath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/nutritionCEST/';
cd(Basepath);
SubjectPath=dir('*KHTCH*');
%%
for i=1:numel(SubjectPath)
    
    
SubjectFullPath=fullfile(SubjectPath(i).folder,SubjectPath(i).name);
cd(SubjectFullPath);
T1MapPath=dir('*Perfusion*');
T1datapath=fullfile(T1MapPath.folder,filesep,T1MapPath.name);
cd(T1datapath);

%pre-load datasets from folder
load('Data.mat');%size(Data)=96x96x1x2x16, 16vtr. 2 for selective and global
load('Map_results.mat');
load('ROIResult_CtxMedulla.mat');
load('AcqAndLsqFit.mat');

Ctx1Mask=roi(1).volume;
Ctx1Mask=rot90(Ctx1Mask,1);
Ctx2Mask=roi(2).volume;
Ctx2Mask=rot90(Ctx2Mask,1);
TI=Data.headerHistory{1}.method.PVM_FairTIR_Arr
%sizes(fitResults)=96x96x1x3x2,1for slice,3 for 3 fitParams, 2 for
% global(1) and selective(2). in this case we would only require global.


for sl=1:size(Data,3)
  for xi=1:size(Data,1)
        for yi=1:size(Data,2)
            if Ctx1Mask(xi, yi)== 0
                continue;
            else
                   FitVector=squeeze(AcqDataGS(xi,yi,sl,:,1));
                   fun=@(A,X)abs(A(1)*(1-2*A(3)*exp(-X./A(2))));%A(3) inversion efficacy
                   Params=fitResults(xi,yi,sl,:,1);%glonal(1),selective(2)
                   plot(TI,FitVector,'*');
                   hold on;
                   plot(linspace(TI(1),TI(end),100),fun(Params,linspace(TI(1),TI(end),100)),'-');
            end

         end
  end
title(sprintf('FAIR T1Map: %s(Left Ctx)',Data.headerHistory{1}.patientID),'FontSize',15)
xlabel('TI time (ms)','FontSize',15);
ylabel('M0 (a.u.)','FontSize',15);
legend({'AcqDataPoint','FitResults'},'FontSize',15,'Location','SouthEast');
grid on;
saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string'))); 
end
close all;




for sl=1:size(Data,3)
  for xi=1:size(Data,1)
        for yi=1:size(Data,2)
            if Ctx2Mask(xi, yi)== 0
                continue;
            else       
                   FitVector=squeeze(AcqDataGS(xi,yi,sl,:,1));
                   fun=@(A,X)abs(A(1)*(1-2*A(3)*exp(-X./A(2))));%A(3) inversion efficacy
                   Params=fitResults(xi,yi,sl,:,1);%glonal(1),selective(2)
                   plot(TI,FitVector,'*');
                   hold on;
                   plot(linspace(TI(1),TI(end),100),fun(Params,linspace(TI(1),TI(end),100)),'-');
            end

        end
   end
 title(sprintf('FAIR T1Map: %s(Right Ctx)',Data.headerHistory{1}.patientID),'FontSize',15)
 xlabel('TI time (ms)','FontSize',15);
 ylabel('M0 (a.u.)','FontSize',15);
 legend({'AcqDataPoint','FitResults'},'FontSize',15,'Location','SouthEast');
 grid on;
 saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
 saveas(gcf, sprintf('%s.fig',get(get(gca,'title'),'string'))); 
end


end





