function [AnimalID,Zref,glucoCESTSpectumValue,Asym_SelectFreq] = WZ_CESTSpectrumInSelectedHz3PtMethod(Datapath,CESTppm,name,metabolite)
%Wan-Ting Zhao 2020-02-12
%this function is customized for renal CEST Asymmetric mesurement of APT* and NOE*
%three point method was first introduced by Kim, S.G. group. Instead of
%conventional MTRasym method, they introduce a reference Z spectrum by
%[Z(ceil)+Z(floor)]./2 to substract Zselected_ppm. Input must include the
%CEST metabolite 'NOE' or 'APT'.

       Data=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
       AnimalID = Data.header.patientID;  

         cd(Datapath);
         CESTROI.name=name;
         if isfile(sprintf('spectrum_%s.mat',CESTROI.name))==1
           load(sprintf('spectrum_%s.mat',CESTROI.name));
         else
            disp('no proper spectrum');
         end

        xaxis4k=linspace(-2000,2000,4000);%
        ROIMeanBase4k=spline(xaxisAcq,ROIMeanAcqBaseCorrected,xaxis4k);
    
%     %debug plot
%     figure,plot(xaxisCorrected,ROIMeanAcqBaseCorrected,'bo',xaxis,ROIMeanBase,'r-')
    
        
        %find minimum as zero and correct
         [min_value,min_idx]=min(ROIMeanBase4k);
         a=xaxis4k(1,min_idx);
         xaxis_corrected1=xaxis4k-a;
            
        pp1=spline(xaxis_corrected1,ROIMeanBase4k);
        splinedROIMeanBase4k=ppval(pp1,xaxis_corrected1);%splineData

         
         GlucoCESTWorkHz = CESTppm*400;
         glucoCESTSpectumValue=ppval(pp1,GlucoCESTWorkHz);
        

        
        if contains(metabolite,'APT')
                GlucoCESTWorkHz = 3.5*400;
                glucoCESTSpectumValue=ppval(pp1,GlucoCESTWorkHz);
                Z4ppm=ppval(pp1,4*400);
                Z3ppm=ppval(pp1,3*400);
                Zref=(Z4ppm+Z3ppm)./2;
                Asym_SelectFreq=Zref-glucoCESTSpectumValue;
        elseif contains(metabolite,'NOE')
                GlucoCESTWorkHz = -3.5*400;
                glucoCESTSpectumValue=ppval(pp1,GlucoCESTWorkHz);
                Zneg2ppm=ppval(pp1,-2*400);
                Z5ppm=ppval(pp1,5*400);%reference said 5ppm instead of -5ppm is better for kidney(Liu, J. DOI:10.1002/nbm.3942)
                Zref=(Zneg2ppm+Z5ppm)./2;
                Asym_SelectFreq=Zref-glucoCESTSpectumValue;                  
        else
                 %calculate Asymetric
                glucoCESTSpectumValue_Positive=ppval(pp1,GlucoCESTWorkHz);
                glucoCESTSpectumValue_Negative=ppval(pp1,-GlucoCESTWorkHz);
                Asym_SelectFreq=glucoCESTSpectumValue_Negative-glucoCESTSpectumValue_Positive;         
        end
        
disp(fprintf('%s AsymCEST %3.2f ppm=%6.4f',AnimalID,CESTppm, Asym_SelectFreq));

cd('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/');
ResultDataPath=fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('AsymCESTppm_%3.2f',CESTppm)));

if isfolder(ResultDataPath)==0

mkdir(sprintf('AsymCESTppm_%3.2f',CESTppm));
else
cd (fullfile('/Users/clairezhao/Documents/MATLAB/CEST/RUNreports/',(sprintf('AsymCESTppm_%3.2f',CESTppm))));
end


ResultDataPathDate=fullfile(ResultDataPath,date);


if isfolder(ResultDataPathDate)==0
    mkdir(ResultDataPath, date);
    cd([ResultDataPath,'/',date]);
%     disp('0')
else
    cd([ResultDataPath,'/',date]);
%     disp('1')
end

save(sprintf('SelectFrqResult_%3.2fppm-%s-%s%s.mat',CESTppm,AnimalID, Data.header.protocolName,Data.headerHistory{1}.acqp.ACQ_time),'GlucoCESTWorkHz', 'Asym_SelectFreq');

end
