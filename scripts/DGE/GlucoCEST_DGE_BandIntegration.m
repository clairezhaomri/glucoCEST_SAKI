function [AnimalID, DeltaS_Map] = GlucoCEST_DGE_BandIntegration(pathBase, pathInfu, ppmRange)
% DGE-style glucoCEST with band integration
% ppmRange: [ppm_min, ppm_max], e.g., [0.5, 2] or [1, 3]
%
% Modified for DGE analysis 2026 Feb

%% Load data (same as version 1)
CESTbase = MRIDataLoader.LoadFrom2DSeq([pathBase filesep 'pdata' filesep '1' filesep '2dseq']);
AnimalID = CESTbase.header.patientID;

if isfile([pathBase filesep 'CESTCorrectData.mat'])
    load([pathBase filesep 'CESTCorrectData.mat'])
    CESTbase(:,:,:) = CESTCorrectData(:,:,:);
end
CESTbase.Permute('echo','read','line');

if strcmp(CESTbase.headerHistory{1}.method.M0Option,'M0_Last')
    S0_base = squeeze(CESTbase.vol(CESTbase.Size('echo'),:,:));
    CESTbase.vol(CESTbase.Size('echo'),:,:) = [];
elseif strcmp(CESTbase.headerHistory{1}.method.M0Option,'M0_First')
    S0_base = squeeze(CESTbase.vol(1,:,:));
    CESTbase.vol(1,:,:) = [];
end

CESTinfu = MRIDataLoader.LoadFrom2DSeq([pathInfu filesep 'pdata' filesep '1' filesep '2dseq']);
if isfile([pathInfu filesep 'CESTCorrectData.mat'])
    load([pathInfu filesep 'CESTCorrectData.mat'])
    CESTinfu(:,:,:) = CESTCorrectData(:,:,:);
end
CESTinfu.Permute('echo','read','line');

if strcmp(CESTinfu.headerHistory{1}.method.M0Option,'M0_Last')
    S0_infu = squeeze(CESTinfu.vol(CESTinfu.Size('echo'),:,:));
    CESTinfu.vol(CESTinfu.Size('echo'),:,:) = [];
elseif strcmp(CESTinfu.headerHistory{1}.method.M0Option,'M0_First')
    S0_infu = squeeze(CESTinfu.vol(1,:,:));
    CESTinfu.vol(1,:,:) = [];
end

MT_Offsets = CESTbase.headerHistory{1}.method.MT_Offsets_NoM0(:).';
localSizeRead = CESTbase.Size('read');
localSizeLine = CESTbase.Size('line');

DeltaS_Map = zeros(localSizeRead, localSizeLine);

%% Define integration range in Hz
freq_min = ppmRange(1) * 400; % Adjust 400 based on your field strength
freq_max = ppmRange(2) * 400;

%% Pixel-wise processing
tic
for k = 1:localSizeRead
    fprintf('Processing pixel %3.3i/%3.3i\n', k, localSizeRead);
    for l = 1:localSizeLine
        
        if S0_base(k,l) == 0, S0_base(k,l) = S0_base(k-1,l); end
        if S0_infu(k,l) == 0, S0_infu(k,l) = S0_infu(k-1,l); end
        
        %% Process baseline
        x_base = MT_Offsets;
        y_base = (CESTbase.vol(:,k,l) ./ S0_base(k,l))';
        [~, xIndex] = sort(x_base);
        CEST_raw_base = [x_base(xIndex); y_base(xIndex)];
        
        [~, min_idx] = min(CEST_raw_base(2,:));
        w0_base = CEST_raw_base(1, min_idx);
        CEST_corrected_base = [CEST_raw_base(1,:) - w0_base; CEST_raw_base(2,:)];
        
        x_spline_base = linspace(CEST_corrected_base(1,1), CEST_corrected_base(1,end), 4000);
        y_spline_base = spline(CEST_corrected_base(1,:), CEST_corrected_base(2,:), x_spline_base);
        [~, min_idx_spline] = min(y_spline_base);
        w0_spline_base = x_spline_base(min_idx_spline);
        x_final_base = x_spline_base - w0_spline_base;
        pp_base = spline(x_final_base, y_spline_base);
        
        %% Process infusion
        x_infu = MT_Offsets;
        y_infu = (CESTinfu.vol(:,k,l) ./ S0_infu(k,l))';
        [~, xIndex] = sort(x_infu);
        CEST_raw_infu = [x_infu(xIndex); y_infu(xIndex)];
        
        [~, min_idx] = min(CEST_raw_infu(2,:));
        w0_infu = CEST_raw_infu(1, min_idx);
        CEST_corrected_infu = [CEST_raw_infu(1,:) - w0_infu; CEST_raw_infu(2,:)];
        
        x_spline_infu = linspace(CEST_corrected_infu(1,1), CEST_corrected_infu(1,end), 4000);
        y_spline_infu = spline(CEST_corrected_infu(1,:), CEST_corrected_infu(2,:), x_spline_infu);
        [~, min_idx_spline] = min(y_spline_infu);
        w0_spline_infu = x_spline_infu(min_idx_spline);
        x_final_infu = x_spline_infu - w0_spline_infu;
        pp_infu = spline(x_final_infu, y_spline_infu);
        
        %% Band integration: trapezoid rule over freq_min to freq_max
        x_integration = linspace(freq_min, freq_max, 100);
        y_base_band = ppval(pp_base, x_integration);
        y_infu_band = ppval(pp_infu, x_integration);
        
        % Mean signal over the band
        S_base_mean = trapz(x_integration, y_base_band) / (freq_max - freq_min);
        S_infu_mean = trapz(x_integration, y_infu_band) / (freq_max - freq_min);
        
        %% Calculate ΔS
        DeltaS_Map(k,l) = (S_base_mean - S_infu_mean) * 100;
        
    end
end
toc

%% Visualization
figure;
imagesc(rot90(DeltaS_Map, 3), [-10 10]);
colormap hot; colorbar;
title(sprintf('DGE Band Map (%.1f-%.1f ppm): %s', ppmRange(1), ppmRange(2), AnimalID));

save([pathBase filesep sprintf('DGE_BandMap_%.1f-%.1fppm.mat', ppmRange(1), ppmRange(2))], 'DeltaS_Map');
saveas(gcf, [pathBase filesep sprintf('DGE_BandMap_%.1f-%.1fppm.fig', ppmRange(1), ppmRange(2))]);

end
