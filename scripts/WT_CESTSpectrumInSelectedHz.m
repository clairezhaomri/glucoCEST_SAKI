function [Asym_SelectFreq] = WT_CESTSpectrumInSelectedHz(Datapath,glucoCESTppm,name)
%this function is designed to spline at selected freqeuncy band and store in  


Asym_SelectFreq = NaN;

try
    Data=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
catch ME
    fprintf('Error: Cannot load data from %s: %s\n', Datapath, ME.message);
    return;
end
        
cd(Datapath);


if isfile(sprintf('spectrum_%s.mat',name))
    try
        load(sprintf('spectrum_%s.mat',name));
        xaxis=linspace(-2000, 2000, 4000);
    catch ME
        fprintf('Error: Cannot load spectrum_%s.mat: %s\n', name, ME.message);
        return;
    end
elseif isfile(sprintf('CESTROIValue_%s.mat',name))
    try
        load(sprintf('CESTROIValue_%s.mat',name));
    catch ME
        fprintf('Error: Cannot load CESTROIValue_%s.mat: %s\n', name, ME.message);
        return;
    end
else
    fprintf('Warning: %s does not have a ROI file (spectrum_%s.mat or CESTROIValue_%s.mat)\n', ...
        Data.header.patientID, name, name);
    return;
end

if ~exist('ROIMeanBase', 'var')
    fprintf('Error: ROIMeanBase variable not found in ROI file for %s\n', Data.header.patientID);
    return;
end

if ~exist('xaxis', 'var')
    fprintf('Error: xaxis variable not found in ROI file for %s\n', Data.header.patientID);
    return;
end

try
    [min_value,min_idx]=min(ROIMeanBase);
    B0Shift=xaxis(1,min_idx);
    xaxis_corrected1=xaxis-B0Shift;
        
    pp1=spline(xaxis_corrected1,ROIMeanBase);
    splinedROIMeanBase4k=ppval(pp1,xaxis);
         
    GlucoCESTWorkHz = glucoCESTppm*400;
    
    glucoCESTSpectumValue_Positive=ppval(pp1,GlucoCESTWorkHz);
    glucoCESTSpectumValue_Negative=ppval(pp1,-GlucoCESTWorkHz);
    
    if glucoCESTppm==3.5 % APT
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
    end
    
    fprintf('Successfully calculated %s glucoCEST %.2f ppm = %.4f\n', ...
        Data.header.patientID, glucoCESTppm, Asym_SelectFreq);
        
catch ME
    fprintf('Error during processing %s at %.2f ppm: %s\n', ...
        Data.header.patientID, glucoCESTppm, ME.message);
    Asym_SelectFreq = NaN;
end

end