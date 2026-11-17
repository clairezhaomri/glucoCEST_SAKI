function [Asym_SelectFreq] = SW_CESTSpectrumInSelectedHz(Datapath,glucoCESTppm,name)
%this function is designed to spline at selected freqeuncy band and store in  
       Data=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
        
        cd(Datapath);
         % ROI
%          CESTROI.name=name;
         
         if isfile(sprintf('spectrum_%s.mat',name))==1
            load(sprintf('spectrum_%s.mat',name));
        
         elseif isfile(sprintf('CESTROIValue_%s.mat',name))==1
            load(sprintf('CESTROIValue_%s.mat',name));
         else
            print(sprintf('%s does not have a ROI',AnimalID));
         end

%         xaxis4k=linspace(-2000,2000,4000);%
% 
% %         ROIMeanBase=pchip(xaxisCorrected,ROIMeanAcqBaseCorrected,xaxis);
%         ROIMeanBase4k=spline(xaxisAcqCorrect,ROIMeanAcqBaseCorrected,xaxis4k);%seems
% %       to generate extra (undesired) data when spline too much, whereas
% %       pchip stays quite fine

    
%     %debug plot
%     figure,plot(xaxisCorrected,ROIMeanAcqBaseCorrected,'bo',xaxis,ROIMeanBase,'r-')
    
        
        %find minimum as zero and correct
        [min_value,min_idx]=min(ROIMeanBase);
         B0Shift=xaxis(1,min_idx);
         xaxis_corrected1=xaxis-B0Shift;
            
        pp1=spline(xaxis_corrected1,ROIMeanBase);
        splinedROIMeanBase4k=ppval(pp1,xaxis_corrected1);%splineData

         
         
        %find the idx
        GlucoCESTWorkHz = glucoCESTppm*400;
        
        %calculate Asymetric
        glucoCESTSpectumValue_Positive=ppval(pp1,GlucoCESTWorkHz);
        glucoCESTSpectumValue_Negative=ppval(pp1,-GlucoCESTWorkHz);
        
        if glucoCESTppm==3.5%APT
            GlucoCESTWorkHz4ppm=4*400;
            GlucoCESTWorkHz3ppm=3*400;

            glucoCESTSpectumValue_4ppm=ppval(pp1,GlucoCESTWorkHz4ppm);
            glucoCESTSpectumValue_3ppm=ppval(pp1,GlucoCESTWorkHz3ppm);
            glucoCESTSpectumValue_APT=ppval(pp1,GlucoCESTWorkHz);
            Asym_SelectFreq=(glucoCESTSpectumValue_4ppm+glucoCESTSpectumValue_3ppm)./2-glucoCESTSpectumValue_APT;
            
        elseif glucoCESTppm==-3.5
            GlucoCESTWorkHzneg2ppm=-2*400;
            GlucoCESTWorkHz5ppm=5*400;
            
            glucoCESTSpectumValue_neg2ppm=ppval(pp1,GlucoCESTWorkHzneg2ppm);
            glucoCESTSpectumValue_5ppm=ppval(pp1,GlucoCESTWorkHz5ppm);
            glucoCESTSpectumValue_NOE=ppval(pp1,GlucoCESTWorkHz);
            Asym_SelectFreq=(glucoCESTSpectumValue_neg2ppm+glucoCESTSpectumValue_5ppm)./2-glucoCESTSpectumValue_NOE;
            
        else
               Asym_SelectFreq=glucoCESTSpectumValue_Negative-glucoCESTSpectumValue_Positive;
    %          Asym_Division=glucoCESTSpectumValue_Negative./glucoCESTSpectumValue_Positive;
        end
        
    
    

% disp(fprintf('glucoCEST %3.2f ppm=%6.4f',glucoCESTppm, Asym_SelectFreq));
% save(sprintf('SelectFrqResult_%s.mat','glucoCESTppm'),'glucoCESTSpectumValue_Positive','glucoCESTSpectumValue_Negative' , 'Asym_SelectFreq','xaxis4k','ROIMeanBase4k');

% saveas(gcf,sprintf('CESTasym:1_2to2_4:%s_%s.png',datetime));
% savefig(gcf, sprintf('CESTasym:1_2to2_4%s.fig',char(datetime))); 
% cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/');
% ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm)));
% 
% if isfolder(ResultDataPath)==0
% 
% mkdir(sprintf('glucoCESTppm_%3.2f',glucoCESTppm));
% else
% cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('glucoCESTppm_%3.2f',glucoCESTppm))));
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
% 
% 
% 
% save(sprintf('SelectFrqResult_%3.2fppm-%s-%s%s.mat',glucoCESTppm,Data.header.patientID, Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'GlucoCESTWorkHz', 'Asym_SelectFreq');

% save(sprintf('spectrum_1_2to2_4%d_%d_%s.mat',StartContrast,EndContrast,datetime),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis');

%additional attribute to daniel folder
% save(sprintf('spectrum_%s_%s_%s.mat',Data.header.patientID,Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'ROIMeanAcqBaseCorrected','xaxisCorrected','xaxis_corrected1','ROIMeanBase','asym_Base','IntegralValue','xaxis');


end
