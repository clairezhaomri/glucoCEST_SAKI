
%% load data
clear all,clc;
% BasePath='/data/projects/bruker-hoerr-cest/data/ANEST-71670__HemeCEST_04_baseline__2019-10-02_20.28.59/';
BasePath='/Users/clairezhao/Documents/data/MRI/pv6/CEST/glucoCEST-coop-heme/ANEST_71671__HemeCEST_05_Heme__2019-10-08_00.56.39/';
cd(BasePath);

%% load ROI to replot
close all,clc;
% load 'ROI_Medulla_R_results.mat'
% load 'ROI_Medulla_L_results.mat'
% load 'ROI_Cortex_L_results.mat'
load 'ROI_Cortex_R_results.mat'

xaxis=linspace(-1800,1800,1000);
  %make symmetric
              [min_value,min_idx]=min(ROIMeanBase);
         a=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-a;
        
        %generate symmetric datapoint
        synData_x_base=linspace(-1800,1800,1000);
        synData_y_base=ROIMeanBase;

        synData_y1_base= synData_y_base;
        FillNr=(length(synData_y_base)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));

        %make sure 500pt in the left hand side
        if FillNr> 0
        synData_y1_base(1:FillNr-1)=synData_y_base(1:FillNr-1);
        synData_y1_base(FillNr:end)=synData_y_base(1:end-FillNr+1);% now symmetric
        synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip

       
        elseif FillNr< 0
        qFillNr=abs(FillNr);       
        synData_y1_base(1:qFillNr)=0;
        synData_y1_base(1:end-qFillNr)=synData_y_base(qFillNr+1:end);
        synData_y1_base(end-qFillNr+1:end)=synData_y_base(end-qFillNr+1:end);% now symmetric
        synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip
        end
        

        %calculate asymSpectrum
        asym_Base=abs(synData_y_base-synData_y1_base);
        
        %find display idx
        start_Contrast = 0;
        startPositions_Contrast1 = find(synData_x_base<start_Contrast);
        startIdx1 = max(startPositions_Contrast1)+1;
%                 
        stop_Contrast = 1800;
        stopPositions_Contrast1 = find(synData_x_base<stop_Contrast);
        stopIdx1 = max(stopPositions_Contrast1)+1;

        %-----------------
        
        min_value=[];
        min_idx=[];
        FillNr=[];
        %-----------------
        [min_value,min_idx]=min(ROIMeanInfu);
         a=xaxis(1,min_idx);
         xaxis_corrected2=xaxis-a; 
        %generate symmetric datapoint
        synData_x_Infu=linspace(-1800,1800,1000);
        synData_y_infu=ROIMeanInfu;

        synData_y1_Infu= synData_y_infu;
        FillNr=(length(synData_y_infu)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));

        %make sure 500pt in the left hand side
        if FillNr> 0
        synData_y1_Infu(1:FillNr-1)=synData_y_infu(1:FillNr-1);
        synData_y1_Infu(FillNr:end)=synData_y_infu(1:end-FillNr+1);% now symmetric
        synData_y1_Infu(end/2+1:end)=flipud(synData_y1_Infu(1:end/2));%now flip

       
        elseif FillNr< 0
        qFillNr=abs(FillNr);       
        synData_y1_Infu(1:qFillNr)=0;
        synData_y1_Infu(1:end-qFillNr)=synData_y1_Infu(qFillNr+1:end);
        synData_y1_Infu(end-qFillNr+1:end)=synData_y1_Infu(end-qFillNr+1:end);% now symmetric
        synData_y1_Infu(end/2+1:end)=flipud(synData_y1_Infu(1:end/2));%now flip
        end
        %calculate asymSpectrum
        asym_Infu=abs(synData_y_infu-synData_y1_Infu);
        %find display idx
        start_Contrast = 0;
        startPositions_Contrast2 = find(synData_x_Infu<start_Contrast);
        startIdx2 = max(startPositions_Contrast2)+1;
%                 
        stop_Contrast = 1800;
        stopPositions_Contrast2 = find(synData_x_Infu<stop_Contrast);
        stopIdx2 = max(stopPositions_Contrast2)+1;
        
        %correct display
%         asym_Base(1:startIdx1)=0;
%         asym_Base(stopIdx1:end)=0;
%         asym_Infu(1:startIdx2)=0;
%         asym_Infu(stopIdx2:end)=0;
  
  figure,errorbar(xaxis_corrected1, ROIMeanBase,ROISemBase,'b');hold on;   
  errorbar(xaxis_corrected2, ROIMeanInfu,ROISemInfu,'r');
  plot(xaxis_corrected1(startIdx1:stopIdx1),asym_Base(startIdx1:stopIdx1), 'b');
  plot(xaxis_corrected2(startIdx2:stopIdx2), asym_Infu(startIdx2:stopIdx2),'r');
%   plot(xaxis_corrected1(startIdx1:stopIdx1),asym_Base(startIdx1:stopIdx1), 'b--');
%   plot(xaxis_corrected2(startIdx2:stopIdx2), asym_Infu(startIdx2:stopIdx2),'r--');
  
  xlabel('Offset (Hz)');
  ylabel('S/S_0');
  axis([-1800 1800 0 1]);
   title('sbj#71671_h_e_m_eCortex_R');
%    title('sbj#71671_b_a_s_eMedulla_R');
   
   
   %% save figure
     saveas(gcf, sprintf('%s.png',get(get(gca,'title'),'string')));
