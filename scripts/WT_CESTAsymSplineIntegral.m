function [AnimalID,IntegralValue,RAW_x,RAW_y] = WT_CESTAsymSplineIntegral(Datapath,StartContrastPpm,EndContrastPpm, name)
%this function is designed to spline at selected freqeuncy band and store in  

% initiation
AnimalID = 'UNKNOWN';
IntegralValue = NaN;
RAW_x = NaN;
RAW_y = NaN;

try
    Data=MRIDataLoader.LoadFrom2DSeq([Datapath filesep 'pdata' filesep '1' filesep '2dseq']); 
    AnimalID = Data.header.patientID;
catch ME
    fprintf('Error: Cannot load data from %s: %s\n', Datapath, ME.message);
    return;
end
          
cd(Datapath);

% check file existance
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
    fprintf('Warning: %s does not have a ROI file (spectrum_%s.mat or CESTROIValue_%s.mat)\n', AnimalID, name, name);
    return;
end

% check variable existance
if ~exist('ROIMeanBase', 'var')
    fprintf('Error: ROIMeanBase variable not found in ROI file for %s\n', AnimalID);
    return;
end

if ~exist('xaxis', 'var')
    fprintf('Error: xaxis variable not found in ROI file for %s\n', AnimalID);
    return;
end

% try asymmetric spectrum
try
    % find minimum for B0 correction
    [min_value,min_idx]=min(ROIMeanBase);
    B0Shift=xaxis(1,min_idx);
    xaxis_corrected1=xaxis-B0Shift;
    
    RAW_x=linspace(xaxis(1),xaxis(end),33);

    % find integral region
    StartContrast1 = StartContrastPpm*400;
    startPositions_Contrast1 = find(xaxis<StartContrast1);
    if isempty(startPositions_Contrast1)
        startIdx = 1;
    else
        startIdx = max(startPositions_Contrast1)+1;
    end
             
    EndContrast1 = EndContrastPpm*400;
    stopPositions_Contrast1 = find(xaxis<EndContrast1);
    if isempty(stopPositions_Contrast1)
        stopIdx = length(xaxis);
    else
        stopIdx = max(stopPositions_Contrast1);
    end
    
    
    pp1=spline(xaxis_corrected1',ROIMeanBase);
    xinterval=linspace(-2000,2000,33);
    RAW_y=ppval(pp1,xinterval);
    
    Z=ppval(pp1,linspace(-2000,2000,4000));

    
    MTRasym=Z(end:-1:1)-Z;
    xbin=zeros(1,4000);
    
   
    startIdx = max(1, min(startIdx, 4000));
    stopIdx = max(1, min(stopIdx, 4000));
    if startIdx > stopIdx
        temp = startIdx;
        startIdx = stopIdx;
        stopIdx = temp;
    end
    
    xbin(startIdx:stopIdx)=1;
    IntegralValue=(trapz(xbin,MTRasym));        
    
    fprintf('Successfully processed %s: Integral glucoCEST %.1f-%.1f ppm = %.4f\n', ...
        AnimalID, StartContrastPpm, EndContrastPpm, IntegralValue);
        
catch ME
    fprintf('Error during processing %s: %s\n', AnimalID, ME.message);
    IntegralValue = NaN;
    RAW_x = NaN;
    RAW_y = NaN;
end

end