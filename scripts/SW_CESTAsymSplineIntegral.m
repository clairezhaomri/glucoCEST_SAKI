function [AnimalID,IntegralValue] = SW_CESTAsymSplineIntegral(Datapath,StartContrastPpm,EndContrastPpm, name)
%this function is designed to spline at selected freqeuncy band and store in  
       Data=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
       
%        ACQ_scan_name = Data.headerHistory{1}.acqp.ACQ_scan_name;
%        ACQ_time = Data.headerHistory{1}.acqp.ACQ_time;
       AnimalID = Data.header.patientID;
              
        cd(Datapath);
        
        
        if isfile(sprintf('spectrum_%s.mat',name))==1
            load(sprintf('spectrum_%s.mat',name));

        
         elseif isfile(sprintf('CESTROIValue_%s.mat',name))==1
            load(sprintf('CESTROIValue_%s.mat',name));

            
         else
            print(sprintf('%s does not have a ROI',AnimalID));
         end
%         xaxis4k=linspace(-2000,2000,4000);%
%         ROIMean4K=spline(xaxis1,ROIMeanAcqBase,xaxis4k);
       
%         RAW_Y=ROIMeanAcqBase;
    
%     %debug plot
%     figure,plot(xaxisAcq,ROIMeanAcqBase,'bo',xaxis,ROIMeanBase,'r-')
    
        %prepare data for aymmetric spectrum
        %make symmetric
       [min_value,min_idx]=min(ROIMeanBase);
         B0Shift=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-B0Shift;
         
       RAW_x=linspace(xaxis(1),xaxis(end),33);
       


        
        %generate symmetric datapoint
%         synData_x_base=linspace(-2000,2000,4000);
%         synData_x_base=linspace(xaxis_corrected1(1),xaxis_corrected1(end),4000);
% 
%         synData_y_base=ROIMeanBase;
% 
%         synData_y1_base= synData_y_base;
%         FillNr=(length(synData_y_base)./2)-min_idx;
% %         FillNr=abs(min_idx-(length(synData_y_base)./2));
%         
% 
%         if isrow(synData_y_base)
%             synData_y_base=synData_y_base';
%         end
%         
%          if isrow(synData_y1_base)
%             synData_y1_base=synData_y1_base';
%         end
% 
%         %make sure half of the total pt in the left hand side
%         if FillNr> 0
%         synData_y1_base(1:FillNr-1)=synData_y_base(1:FillNr-1);
%         synData_y1_base(FillNr:end)=synData_y_base(1:end-FillNr+1);% now symmetric
%         synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip
% 
%        
%         elseif FillNr< 0
%         qFillNr=abs(FillNr);       
%         synData_y1_base(1:qFillNr)=0;
%         synData_y1_base(1:end-qFillNr)=synData_y_base(qFillNr+1:end);
%         synData_y1_base(end-qFillNr+1:end)=synData_y_base(end-qFillNr+1:end);% now symmetric
%         synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip
%         
%         elseif FillNr == 0
%         synData_y1_base(1:end/2)=synData_y_base(1:end/2);
%         synData_y1_base(end/2+1:end)=flipud(synData_y_base(1:end/2));
%         end
%         


        %find integral interval of *mirror-ed* data
        StartContrast1 = StartContrastPpm*400;
        startPositions_Contrast1 = find(xaxis<StartContrast1);
        startIdx = max(startPositions_Contrast1)+1;
                 
        EndContrast1 = EndContrastPpm*400;
        stopPositions_Contrast1 = find(xaxis<EndContrast1);
        stopIdx = max(stopPositions_Contrast1);
        
%         %debug plot
%         figure,plot(xaxis_corrected1,synData_y_base,'r');hold on;%synthesized data
%         plot(xaxis,synData_y1_base,'b')%mirror data
%         
        
        pp1=spline(xaxis_corrected1',ROIMeanBase);
%         pp2=spline(xaxis', synData_y1_base);

        xinterval=linspace(-2000,2000,33);
%         y1=ppval(pp1,xaxis);
%         y2=ppval(pp2,xaxis);
        
        RAW_y=ppval(pp1,xinterval);
%         %debug plot        
%         figure, plot(xaxis,y1,'r',xaxis,y2,'b')
        Z=ppval(pp1,linspace(-2000,2000,4000));

        %calculate asymSpectrum
        MTRasym=Z(end:-1:1)-Z;
%         asym_Base=(y2-y1);%mirror trace-signal trace     
%         IntegralValue=abs(trapz(asym_Base(startIdx:stopIdx)));
        xbin=zeros(1,4000);
        xbin(startIdx:stopIdx)=1;
        IntegralValue=(trapz(xbin,MTRasym));        
        %debug plot        
%         figure, plot(xaxis,y1,'r',xaxis,y2,'b')
%         plot(xaxis,asym_Base, 'm-','Linewidth',2);
        
    
    
%     disp(fprintf('Integral glucoCEST %3.1f-%3.1f ppm=%6.4f',StartContrastPpm,EndContrastPpm, IntegralValue))
%     save(sprintf('SelectFrqResult_%s.mat','glucoCESTppm'),'StartContrast','EndContrast' , 'IntegralValue','xaxis1','ROIMeanBase4k');
    

  h=figure; 
  set(h,'Visible','off');%dont show up when queuing
  axes1 = axes('Parent',h);
  
 
  plot(RAW_x./400,RAW_y,'b*','Linewidth',2); hold on;
  plot((xaxis_corrected1)./400, ROIMeanBase,'b-','Linewidth',2)
  ylabel('S/S_0','FontSize',15,'FontWeight','bold');
  xlabel('Offset (ppm)','FontSize',15,'FontWeight','bold');hold on;
  
  
%   if stopIdx>4000
 plot(xaxis(end/2+1:end)./400,MTRasym(end/2+1:end), 'm-','Linewidth',2);   
%   else
%   plot((synData_x_base(startIdx:stopIdx))./400,asym_Base, 'm-','Linewidth',2);
%   end
      
%   title(sprintf('ROI-based CEST & Asymmetric Spectrum(%s)',ROI.title),'FontSize',15);
  title(sprintf('CESTMTRAsym:%s',AnimalID),'FontSize',15);
  legend({'AcqData','Spline','AsymSpectrum'},'FontSize',12, 'Location', 'northeast');
  legend('boxoff')
  
  box(axes1,'on');
  grid on;


    axis([-5 5 0 1]);
    set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [-5 -4 -3 -2 -1 0 1 2 3 4 5],'YTick',...
    [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],'ZTick',[-1 0 1]);
    set(gca, 'XDir','reverse');
    axes('Position',[.66 .15 .2 .2])
    box on
    plot(xaxis./400,MTRasym.*100,'b-','Linewidth',2)
    hold on;
    plot([0,5],[0,0], 'k-')
    axis([0 2.5 -8 8])
    set(gca, 'XDir','reverse')
    ylabel('MTRasym(%)')
    grid on
    
 saveas(gcf,sprintf('WT-CESTasym:%3.1f-%3.1f_%s.fig',StartContrastPpm,EndContrastPpm, AnimalID));


% disp(fprintf('IntegralValue=%6.4f',IntegralValue))

% saveas(gcf,sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.png',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name, ACQ_time));
% savefig(gcf, sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.fig',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name,ACQ_time)); 
% save(sprintf('spectrum:%3.1f-%3.1f_%s.mat',StartContrastPpm,EndContrastPpm,AnimalID,ACQ_scan_name, ACQ_time),'ROIMeanAcqBase','xaxisAcqCorrect','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis1');



% 
% cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/');
% ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm)));
% 
% if isfolder(ResultDataPath)==0
% mkdir(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm));
% cd(ResultDataPath)
% %     disp('0')
% else
% cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.1f_%3.1f',StartContrastPpm,EndContrastPpm))));
% %     disp('1')
% end
% 
% 
% ResultDataPathDate=fullfile(ResultDataPath,date);
% 
% 
% if isfolder(ResultDataPathDate)==0
%     mkdir(ResultDataPath, date);
%     cd([ResultDataPath,'/',date]);
% %     disp('0')
% else
%     cd([ResultDataPath,'/',date]);
% %     disp('1')
% end


% saveas(gcf,sprintf('CESTasym:%3.1f-%3.1f_%s.png',StartContrastPpm,EndContrastPpm, AnimalID));
% savefig(gcf, sprintf('CESTasym:%3.1f-%3.1f_%s%s%s.fig',StartContrastPpm,EndContrastPpm, AnimalID,ACQ_scan_name,ACQ_time)); 
% save(sprintf('spectrum:%3.1f-%3.1f_%s%s%s.mat',StartContrastPpm,EndContrastPpm,AnimalID,ACQ_scan_name, ACQ_time),'ROIMeanAcqBase','xaxisAcqCorrect','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis1');
% close all;
% %additional attribute to daniel folder
% cd ('//data/projects/bruker-hoerr-cest/documents/spline0.5-2ppm/');
% save(sprintf('spectrum_%s_%s_%s.mat',Data.header.patientID,Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'ROIMeanAcqBase','xaxisAcqCorrect','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis');
end
