% -------------------------------------------------------------------------

% MATLAB script adapting MZaiss(2014, NMRbiomed) open source matlab script
% for diet+S-AKI project, Wan-Ting Zhao Ver.20250902
% -------------------------------------------------------------------------

%% load data
% Data = MRIDataLoader.LoadFrom2DSeq('path/to/2dseq');

BasePath ='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/20210713_070810_KHTCH_11038_1_1CtrlProteinFastAcute/CEST_base/';
Data = MRIDataLoader.LoadFrom2DSeq('/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/20210713_070810_KHTCH_11038_1_1CtrlProteinFastAcute/CEST_base/pdata/1');

Data = Data.Permute('read,line,slice,echo');
% Data = Data.Squeeze();
Data_raw = Data(:,:,:,1:end-1);
M0_stack = Data(:,:,:,end);

%% input data
addpath(genpath(cd('.')))

% we assume that you have the normalization stack M0
% the saturated  stack Mz
% and the parameter struct P containing
P.SEQ.w= linspace(-5, 5, 33);%start offets, stop offsets, number of offsets(exclude S0).
P.EVAL.w_interp=P.SEQ.w;
%%==added 20250820 Wan-Ting Zhao,reordering Mz_stack for Bruker Data===
% load raw bruker data,manually create Data_raw(which is 1:end-1 offsets) and M0_stack.
% Data = MRIDataLoader.LoadFrom2DSeq('path/to/2dseq');
% Data = Data.Permute('read,line,slice,echo');
% Data_raw = Data(:,:,:,1:end-1);
% M0_stack = Data(:,:,:,end);

sorted_w = P.SEQ.w; 
num_offsets = length(sorted_w);
original_w = zeros(1, num_offsets);
neg_offsets = sorted_w(sorted_w < 0);     
pos_offsets = sorted_w(sorted_w >= 0); 

for i = 1:num_offsets/2
    original_w(2*i - 1) = neg_offsets(i);               % e.g., index 1, 3, 5...
    original_w(2*i)     = pos_offsets(end - i + 1);     % e.g., index 2, 4, 6...
end
[~, reorder_indices] = ismember(sorted_w, original_w);
Data_reordered = Data_raw(:,:,:,reorder_indices);
Mz_stack =Data_reordered;
%======================================================================



P.SEQ.stack_dim=size(Mz_stack);

if ~exist('M0_stack') || ~exist('Mz_stack') || ~exist('P')|| ~isfield(P.SEQ,'w')|| ~isfield(P.EVAL,'w_interp')
    errordlg('missing data');
end

%% DEFINE 2D Segment of ones and NaNs
% this creates the mask Segment
Segment= M0_stack>0.9*mean(M0_stack(:,:,1)); %empirical weighting. 0.2 for phantom, 0.8 for sepsis kidney project
% Segment= M0_stack>100;%for phantom intensity thresholding

%% calculate internal dB0 map
[dB0_stack_int yS] = MINFIND_SPLINE_3D(Mz_stack,Segment,P);

%% correct data for dB0 with internal B0-map from the measurement
[Mz_CORR ] = B0_CORRECTION(Mz_stack,dB0_stack_int,P,Segment);

[Z_uncorr] = NORM_ZSTACK(Mz_stack,M0_stack,P,Segment);
%% Normalize uncorrected data and B0 corrected data
[Z_corrExt] = NORM_ZSTACK(Mz_CORR,M0_stack,P,Segment);


%% start imgui 
% close(imgui); imgui %only MTRasym spline at this point

%% Multi-Lorentzian fitting

P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];%WTZ: options of lsqcurvefit, e.g.tolerance of functional value, tolerance of step size,when to stop the algorithm when it finds a local minimum..etc.
P.FIT.nIter     = 50;
P.FIT.modelnum  = 6;     %% number of Lorentzian pools

P.FIT.extopt=1; % change parameters explicitly
%Lorentzian line LI defined by amplitude Ai, width Gi [ppm] and offset dwi[ppm]: Li=Ai.*Gi^2/4./ (Gi^2/4+(dw-dwi).^2) ;
%1=water; 2=amide; 3=NOE; 4=MT; 5=amine; 6=glucose
%const.Zi   A1    G1    dw1     A2     G2    dw2     A3    G3   dw3    A4
%G4   dw4    A5    G5   dw5     A6     G6    dw6
lb = [ 0.5  0.02  0.3  -1       0.0    0.4   +3     0.0    1    -4.5    0.0   10   -4     0.0   1    1   0.0   0.1   0.0 ];%WTZ: added glucose lb (DOI: https://doi.org/10.14316/pmp.2017.28.4.135)
ub = [ 1    1     10   +1       0.2    3     +4     0.4    5    -2        1   100  -1     0.2   3.5  2.5 1.0   1.8   1.8];%WTZ: added glucose ub
p0 = [ 1    0.9   1.4   0       0.025  0.5   3.5     0.02   3    -3.5    0.1   25   -1     0.01  1.5  2.2 0.02  0.2   0.9];%WTZ: added initial value of glucose pool
P.FIT.lower_limit_fit = lb; 
P.FIT.upper_limit_fit = ub; 
P.FIT.start_fit = p0;

%% fitting
tic
[popt, P] = FIT_3D(Z_corrExt,P,Segment); % perform the fit pixelwise
toc
%%
[Zlab, Zref] = get_FIT_LABREF(popt,P,Segment);  % create Reference values Z_Ref=(Z_lab - Li)

%   imgui


%% debug glucose map
% glucoseAmp = squeeze(popt(:,:,1,17));
% figure,imagesc(fliplr(flipud(rot90(glucoseAmp .* Segment(:,:,1),1))),[0,0.15]),colormap hot, colorbar('FontSize',18),axis image, axis off
% title('Infusion Amplitude of Glc.','FontSize',18)
%% debug glucose pool integral map
glucose_AmpMap = squeeze(popt(:,:,1,17));
glucose_GMap = squeeze(popt(:,:,1,18));
gluc_integral_fullMap = glucose_AmpMap .* glucose_GMap .* pi/2;   % units: amplitude × ppm, integral the entire glc pool

figure,imagesc(fliplr(gluc_integral_fullMap .* Segment(:,:,1)),[0,0.2]),colormap jet, colorbar('FontSize',18),axis image, axis off %phantom
title('Integral of full glucose pool','FontSize',18)

% mask and result
test =gluc_integral_fullMap.*ROI;
figure,imagesc(test,[0,0.2]),colormap jet

%% ROI-based analysis
load('ROI_ctx.mat');% at BasePath level.
ROI_rot = rot90(ROI, 3);%rotation for visualize
ROI_nan_rot = double(ROI_rot);
ROI_nan_rot(ROI_nan_rot == 0) = NaN;

M0_rot = rot90(M0_stack, 3);
gluc_rot = rot90(gluc_integral_fullMap, 3);
croppedMap_rot = gluc_rot .* ROI_nan_rot;

% Figure 1 — mask
figure('Position', [100, 100, 500, 500])
imagesc(M0_rot, [0, 150])
colormap gray
axis image
hold on

roi_red = cat(3, ones(size(ROI_rot)), zeros(size(ROI_rot)), zeros(size(ROI_rot)));
h_roi = imshow(roi_red);
set(h_roi, 'AlphaData', 0.3 * double(ROI_rot > 0));
title('mask', 'FontSize', 18);
hold off

set(gca, 'Position', [0.05, 0.05, 0.75, 0.9]);  % SAME box as figure 2

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, 'mask.png', '-dpng', '-r300');


% Figure 2 - Glucose overlay
figure('Position', [100, 100, 500, 500])

% --- Background axes ---
ax1 = axes;
axes(ax1);
imagesc(M0_rot, [0, 150])
colormap(ax1, gray)
axis image
axis off

% --- Overlay axes ---
ax2 = axes;
axes(ax2);
h_roi2 = imagesc(croppedMap_rot);
colormap(ax2, jet)
axis image
set(h_roi2, 'AlphaData', 0.7 * double(~isnan(croppedMap_rot)));

caxis(ax2, [0, 0.2]);   % pre-R2022a equivalent of clim

cb = colorbar(ax2);
cb.FontSize = 18;
cb.Label.String = 'Glucose integral';

ax2.XColor = 'none';
ax2.YColor = 'none';
ax2.Color = 'none';

linkaxes([ax1, ax2])
ax1.Position = ax2.Position;

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, 'glucose_roi_overlay.png', '-dpng', '-r300');


% Compute mean and std excluding NaN
roi_mean = mean(croppedMap_rot(:), 'omitnan');
roi_std  = std(croppedMap_rot(:), 'omitnan');

fprintf('Mean: %.4f, SD: %.4f\n', roi_mean, roi_std);


%% debug fitted lorentzian multi-pool curve include glucose
% Step 1: Extract mean Z-spectrum from ROI
N_offsets = size(Z_corrExt, 4);
meanZspec = zeros(N_offsets, 1);

for jj = 1:N_offsets
    slice_img = squeeze(Z_corrExt(:,:,1,jj));
    voxels    = slice_img(ROI == 1);
    meanZspec(jj) = mean(voxels(~isnan(voxels)));
end

% Step 2: Fit the mean Z-spectrum
[P] = fitmodelfunc_NUM(meanZspec, P);
popt_roi = levmar_fit(meanZspec, P);

% Step 3: Reconstruct all pool functions
[f, fZi, f1, f2, f3, f4, f5, f6] = fitmodelfunc_ANA(popt_roi, P);

% Step 4: Plot with yyaxis
w_plot = linspace(min(P.SEQ.w), max(P.SEQ.w), 500)';
W      = P.SEQ.w;

figure('Color','white','Position',[100 100 900 550]);

% --- LEFT Y-AXIS: all pools except glucose ---
yyaxis left;
hold on;

h1 = plot(W,        meanZspec,    'ko',  'MarkerFaceColor','k', 'MarkerSize',5);
h2 = plot(w_plot,   f(w_plot),    '-k',  'LineWidth',2);
h3 = plot(w_plot,   f1(w_plot),   '--c', 'LineWidth',1.5);
h4 = plot(w_plot,   f2(w_plot),   '--b', 'LineWidth',1.5);
h5 = plot(w_plot,   f3(w_plot),   '--',  'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
h6 = plot(w_plot,   f4(w_plot),   '--g', 'LineWidth',1.5);
h7 = plot(w_plot,   f5(w_plot),   '--r', 'LineWidth',1.5);

ylabel('Z / Z_0', 'FontSize',12, 'Color','k');
set(gca, 'YColor','k');
ylim([0 1]);

% --- RIGHT Y-AXIS: glucose pool only ---
yyaxis right;

h8 = plot(w_plot, f6(w_plot), '-m', 'LineWidth',2.5);

ylabel('Glucose pool amplitude', 'FontSize',25, 'Color','m');
set(gca, 'YColor','m');
% scale right axis to match glucose signal range (adjust if needed)
ylim([0 0.2]);   % e.g. 0-5%, adjust based on your data

hold off;

set(gca, 'FontSize', 25,'XDir','reverse');
xlim([-4 4]);
xlabel('Offset (ppm)', 'FontSize',25);
% title('Mean ROI Z-spectrum with Lorentzian Fits', 'FontSize',25);
grid on;
box on;

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'MeanROI_Zspectrum_with_LorentzFit.png', '-dpng', '-r300');

%%
close all;
save('metadata'); 
%% Extract Amplitude and LW of the targeted metabolites
% T = extract_pool_perROI(popt, ROI, {'water','amide','NOE','MT','amine','glucose'});
%% define mock T1map if no T1map in identical geometry
% 
% mapSize = 128; 
% 
% % T1 values are converted from ms to seconds.
% T1_cortex = 1880 / 1000;  % 1.793 s
% T1_medulla = 2480 / 1000; % 2.082 s
% T1_background = 1.5;      % A default value for surrounding tissue
% 
% % Initialize the map with the background T1 value
% T1map = ones(mapSize, mapSize) * T1_background;%uniform backgroup
% T1map = T1map.*Segment;%intensity thresholding mask
% 
% % has to be at data directory
% try
%     load('ROI_ctx.mat', 'ROI');
%     cortexMask = logical(ROI); % Ensure the mask is logical (true/false)
%     T1map(cortexMask) = T1_cortex;
% catch ME
%     warning('Could not load or apply ROI_ctx.mat. A mock mask will not be created. Error: %s', ME.message);
% end
% 
% try
%     load('ROI_medulla.mat', 'ROI');
%     medullaMask = logical(ROI); % Ensure the mask is logical
%     T1map(medullaMask) = T1_medulla;
% catch ME
%     warning('Could not load or apply ROI_medulla.mat. A mock mask will not be created. Error: %s', ME.message);
% end
% 
% % mock T1map done.
% 
% %% -------------------------------------------------------------------------
% % pH Map Calculation from an AREX Map
% % -------------------------------------------------------------------------
% % This script first calculates an Apparent Exchange-dependent Relaxation (AREX)
% % map from Z-spectrum and T1 data, and then converts it to a pH map using
% % the formula from Zaiss et al. (2014).
% %
% % This version is configured for 128x128 maps and uses loaded ROI masks
% % to generate a T1 map with specific renal cortex and medulla values.
% %
% % AREX Formula:
% % AREX = (1/Zlab - 1/Zref) / T1map
% %
% % pH Formula (Equation 19):
% % pH = 6.4 + log10( (AREX/DC) / (fb * 5.57) )
% %
% % -------------------------------------------------------------------------
% 
% 
% % --- Step 1: Define Constants based on the paper ---
% 
% % For continuous wave (CW) CEST, the duty cycle (DC) is 1.
% DC = 1;
% 
% % Proton fraction (fb) for in vivo Amide Proton Transfer (APT).
% % This value is reported in the paper for APT in the brain.
% % IMPORTANT: You may need to change this value depending on the tissue
% % and the specific CEST agent you are imaging.
% fb = 0.00115; % Corresponds to 0.115%
% K_CONSTANT = 5.57;
% 
% 
% % focus on amide
% 
% % AREX_map_amide=(1./Zlab - 1./Zref.Amide) ./ T1map.vol;
% AREX_map_amide=(1./Zlab - 1./Zref.Amide) ./ T1map;
% 
% ratio_term = (AREX_map_amide / DC) / (fb * K_CONSTANT);
% pH_maps = nan(size(AREX_map_amide));
% positive_indices = ratio_term > 0;
% pH_maps(positive_indices) = 6.4 + log10(ratio_term(positive_indices));
% 
% 
% %% plot AREX and pH map
% % Display the calculated AREX map
% % pH_maps_amide=squeeze(pH_maps(:,:,1,28));
% % pH_maps_amide=squeeze(pH_maps(:,:,5,40));
% pH_maps_amide=squeeze(pH_maps(:,:,1,35));
% 
% 
% 
% 
% figure,
% % imagesc(rot90(squeeze(AREX_map_amide(:,:,1,28)),3));
% imagesc(rot90(squeeze(AREX_map_amide(:,:,1,35)),3));
% 
% % title('AREX','FontSize',16);
% axis image;
% axis off;
% cb =colorbar;
% set(cb,'FontSize',25);
% ylabel(cb,'AREX');
% colormap(gca, 'magma');
% % caxis([0 0.15]); % Set color limits for better contrast
% caxis([0 0.5]); % Set color limits for better contrast
% 
% % saveas(gca,'AREX-Map.fig');
% % saveas(gca,'AREX-Map.png');
% 
% % Display the resulting pH map
% figure,
% imagesc(rot90(pH_maps_amide,3));
% % title('pH Map','FontSize',16);
% axis image;
% axis off;
% cb =colorbar;
% set(cb,'FontSize',25);
% ylabel(cb,'pH');
% colormap(gca, 'parula');
% % caxis([6.5 8]); % Adjust the color axis for better pH visualization
% caxis([4.4 4.8]); % Adjust the color axis for better pH visualization
% 
% saveas(gca,'pH-Map.fig');
% saveas(gca,'pH-Map.png');
% 
% 
% 
% %% prepare 1d cell for violin plot(function need to download from online)
% pH_1d_ctx = pH_maps_amide(cortexMask);
% pH_1d_medulla = pH_maps_amide(medullaMask);
% % all_ph_values =[pH_1d_ctx,pH_1d_medulla];
% % labels_ctx = repmat({'ctx'}, length(pH_1d_ctx), 1);
% % labels_medulla = repmat({'medulla'}, length(pH_1d_medulla), 1);
% % all_labels = [labels_ctx; labels_medulla];
% 
% figure('Position', [150, 150, 600, 500]);
% violin({pH_1d_ctx, pH_1d_medulla}, ...
%        'xlabel', {'Cortex', 'Medulla'}, ... 
%        'ShowData', true, ...
%        'ShowMedian', true);
%    
% % title('Distribution of pH Values in Different ROIs', 'FontSize', 16);
% ylim([4.5,9]);
% set(gca, 'FontSize', 25);
% % xlabel('Region of Interest', 'FontSize', 16);
% ylabel('pH', 'FontSize', 25);
% legend off;
% grid on;
% saveas(gca,'violin.fig');
% saveas(gca,'violin.png');
% 
% %% 
% % =========================================================================
% % Display Statistics
% % =========================================================================
% 
% % --- Calculate and Display Statistics ---
% % Calculate standard deviation for each region
% mean_ctx = mean(pH_1d_ctx);
% mean_medulla = mean(pH_1d_medulla);
% std_ctx = std(pH_1d_ctx);
% std_medulla = std(pH_1d_medulla);
% median_ctx = median(pH_1d_ctx);
% median_medulla = median(pH_1d_medulla);
% 
% % Display the results in the command window
% fprintf('\n--- Statistical Analysis ---\n');
% fprintf('Cortex pH Mean: %.4f\n', mean_ctx);
% fprintf('Cortex pH Median: %.4f\n', median_ctx);
% fprintf('Cortex pH Std Dev:  %.4f\n', std_ctx);
% fprintf('----------------------------\n');
% fprintf('Medulla pH mean: %.4f\n', mean_medulla);
% fprintf('Medulla pH Median: %.4f\n', median_medulla);
% fprintf('Medulla pH Std Dev:  %.4f\n', std_medulla);
% fprintf('----------------------------\n');
