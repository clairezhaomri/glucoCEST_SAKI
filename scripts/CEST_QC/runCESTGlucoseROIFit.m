function [popt, P, gluc_integral_fullMap, M0_stack, Segment, Z_corrExt] = runCESTGlucoseROIFit(dataDir, varargin)
% RUNCESTGLUCOSEROIFIT  Runs the expensive CEST pixelwise fit ONLY, and
%   saves the results needed for later plotting (no figures generated here).
%
%   [popt, P, gluc_integral_fullMap, M0_stack, Segment, Z_corrExt] = ...
%       runCESTGlucoseROIFit(dataDir)
%
%   Saves LorentzianFitResult.mat into dataDir containing: popt, P,
%   gluc_integral_fullMap, M0_stack, Segment, Z_corrExt.
%
%   NOTE: this function ALWAYS performs the fit when called -- the
%   "skip if LorentzianFitResult.mat already exists" logic lives in the
%   calling queue script (WZ_QueueCEST_Analyze.m), not here.
%
%   CHANGE FROM PREVIOUS VERSION: before loading raw 2dseq data, this now
%   checks dataDir (the CEST_base/CEST_infu level, NOT the pdata
%   subfolder) for a 'CESTCorrectData.mat' file. If present, its
%   'CESTCorrectData' variable is used INSTEAD of the raw
%   MRIDataLoader.LoadFrom2DSeq load -- e.g. for subjects where the raw
%   2dseq data needed manual motion/artifact correction saved separately.
%   If absent, behavior is unchanged (loads raw 2dseq as before).
%
%   IMPORTANT: the number of saturation offsets is still auto-detected
%   from the loaded data itself (size(Data_raw, 4)) rather than assumed
%   fixed at 33 -- this applies whether Data came from the raw 2dseq load
%   or from CESTCorrectData.mat, since the offset count is read AFTER
%   this substitution happens.
%
%   NAME-VALUE PAIRS: same fit-related options as before
%   'PdataNum', 'Woffsets', 'WoffsetsRange', 'SegmentFactor', 'ModelNum',
%   'FitOptions', 'FitNIter', 'Lb', 'Ub', 'P0', 'AddPathDataDir'
%
%     'Woffsets'      - offset vector; if empty (default) or length
%                        mismatches the data, it's auto-generated from
%                        'WoffsetsRange' and the actual detected offset count.
%     'WoffsetsRange' - [min, max] ppm range used when auto-generating
%                        the offset vector (default [-5, 5]).
warning('off', 'curvefit:prepareFittingData:sizeMismatchInVariables');
p = inputParser;
addParameter(p, 'PdataNum', 1);
addParameter(p, 'Woffsets', []);
addParameter(p, 'WoffsetsRange', [-5, 5]);
addParameter(p, 'SegmentFactor', 0.8);
addParameter(p, 'ModelNum', 6);
addParameter(p, 'FitOptions', [1E-04, 1E-15, 1E-10, 1E-04, 1E-06]);
addParameter(p, 'FitNIter', 50);
addParameter(p, 'Lb', [ 0.5 0.02 0.3 -1 0.0 0.4 +3 0.0 1 -4.5 0.0 10 -4 0.0 1 1 0.0 0.5 0.9 ]);
addParameter(p, 'Ub', [ 1 1 10 +1 0.2 3 +4 0.4 5 -2 1 100 -1 0.2 3.5 2.5 1.0 1.8 1.5]);
addParameter(p, 'P0', [ 1 0.9 1.4 0 0.025 0.5 3.5 0.02 3 -3.5 0.1 25 -1 0.01 1.5 2.2 0.05 0.5 1.2]);
addParameter(p, 'AddPathDataDir', false);
parse(p, varargin{:});
opt = p.Results;

if ~isfolder(dataDir)
    error('runCESTGlucoseROIFit:badDir', 'Directory does not exist: %s', dataDir);
end
if opt.AddPathDataDir
    addpath(genpath(dataDir));
end

%% -------- load data --------
pdataPath = fullfile(dataDir, 'pdata', num2str(opt.PdataNum));
if ~isfolder(pdataPath)
    error('runCESTGlucoseROIFit:noPdata', 'pdata folder not found: %s', pdataPath);
end

Data = MRIDataLoader.LoadFrom2DSeq(pdataPath);

% If a corrected data file exists at the CEST_base/CEST_infu level (NOT
% inside pdata), use it INSTEAD of the raw 2dseq load. Otherwise Data
% stays exactly as loaded above.
correctDataPath = fullfile(dataDir, 'CESTCorrectData.mat');
if isfile(correctDataPath)
    fprintf('Found CESTCorrectData.mat in %s -- using corrected data instead of raw 2dseq load.\n', dataDir);
    S_corr = load(correctDataPath, 'CESTCorrectData');
    if ~isfield(S_corr, 'CESTCorrectData')
        error('runCESTGlucoseROIFit:missingCorrectDataVar', ...
            '%s exists but does not contain a variable named CESTCorrectData.', correctDataPath);
    end
    Data = S_corr.CESTCorrectData;
end

Data = Data.Permute('read,line,slice,echo');
Data_raw = Data(:,:,:,1:end-1);
M0_stack = Data(:,:,:,end);

%% -------- auto-detect actual offset count from the data --------
N_offsets_actual = size(Data_raw, 4);

if isempty(opt.Woffsets)
    sorted_w = linspace(opt.WoffsetsRange(1), opt.WoffsetsRange(2), N_offsets_actual);
    fprintf('Woffsets not specified: auto-generated %d-point offset vector from %.1f to %.1f ppm.\n', ...
        N_offsets_actual, opt.WoffsetsRange(1), opt.WoffsetsRange(2));
elseif numel(opt.Woffsets) ~= N_offsets_actual
    warning('runCESTGlucoseROIFit:offsetMismatch', ...
        ['Specified Woffsets has %d points but the loaded data has %d offsets. ' ...
         'Regenerating a %d-point offset vector over [%.1f, %.1f] ppm to match the data. ' ...
         'Verify this matches your actual acquisition protocol for %s.'], ...
        numel(opt.Woffsets), N_offsets_actual, N_offsets_actual, ...
        opt.WoffsetsRange(1), opt.WoffsetsRange(2), dataDir);
    sorted_w = linspace(opt.WoffsetsRange(1), opt.WoffsetsRange(2), N_offsets_actual);
else
    sorted_w = opt.Woffsets;
end

%% -------- parameter struct / offsets --------
P = struct();
P.SEQ.w = sorted_w;
P.EVAL.w_interp = P.SEQ.w;

num_offsets = numel(sorted_w);
neg_offsets = sorted_w(sorted_w < 0);
pos_offsets = sorted_w(sorted_w >= 0);
nNeg = numel(neg_offsets);
nPos = numel(pos_offsets);

% Generalized interleave: works for even OR odd offset counts, and for
% unequal neg/pos counts (leftover elements appended at the end rather
% than crashing on non-integer indexing).
original_w = zeros(1, num_offsets);
idx = 1;
maxPairs = min(nNeg, nPos);
for i = 1:maxPairs
    original_w(idx) = neg_offsets(i); idx = idx + 1;
    original_w(idx) = pos_offsets(end - i + 1); idx = idx + 1;
end
if nNeg > maxPairs
    remaining = neg_offsets(maxPairs+1:end);
    original_w(idx:idx+numel(remaining)-1) = remaining;
    idx = idx + numel(remaining);
end
if nPos > maxPairs
    remaining = pos_offsets(1:end-maxPairs);
    original_w(idx:idx+numel(remaining)-1) = remaining;
end

[~, reorder_indices] = ismember(sorted_w, original_w);
Data_reordered = Data_raw(:,:,:,reorder_indices);
Mz_stack = Data_reordered;

P.SEQ.stack_dim = size(Mz_stack);

if ~exist('M0_stack','var') || ~exist('Mz_stack','var') || ~exist('P','var') ...
        || ~isfield(P.SEQ,'w') || ~isfield(P.EVAL,'w_interp')
    error('runCESTGlucoseROIFit:missingData', 'missing data');
end

% Sanity check: reordering must not have dropped/duplicated offsets
if any(reorder_indices == 0) || numel(unique(reorder_indices)) ~= num_offsets
    error('runCESTGlucoseROIFit:reorderFailed', ...
        'Offset reordering failed for %s (%d offsets). Check acquisition protocol.', ...
        dataDir, num_offsets);
end

%% -------- segmentation mask --------
Segment = M0_stack > opt.SegmentFactor * mean(M0_stack(:,:,1));

%% -------- B0 correction and normalization --------
[dB0_stack_int, ~] = MINFIND_SPLINE_3D(Mz_stack, Segment, P);
[Mz_CORR] = B0_CORRECTION(Mz_stack, dB0_stack_int, P, Segment);
[Z_corrExt] = NORM_ZSTACK(Mz_CORR, M0_stack, P, Segment);

%% -------- multi-Lorentzian fit setup --------
P.FIT.options = opt.FitOptions;
P.FIT.nIter = opt.FitNIter;
P.FIT.modelnum = opt.ModelNum;
P.FIT.extopt = 1;
P.FIT.lower_limit_fit = opt.Lb;
P.FIT.upper_limit_fit = opt.Ub;
P.FIT.start_fit = opt.P0;

%% -------- pixelwise fit (the expensive step) --------
tic
[popt, P] = FIT_3D(Z_corrExt, P, Segment);
toc

%% -------- glucose pool integral map --------
glucose_AmpMap = squeeze(popt(:,:,1,17));
glucose_GMap = squeeze(popt(:,:,1,18));
gluc_integral_fullMap = glucose_AmpMap .* glucose_GMap .* pi/2;

%% -------- save everything needed for plotting later --------
save(fullfile(dataDir, 'LorentzianFitResult.mat'), ...
    'popt', 'P', 'gluc_integral_fullMap', 'M0_stack', 'Segment', 'Z_corrExt');

fprintf('Fit complete (%d offsets). Results saved to %s\n', num_offsets, fullfile(dataDir, 'LorentzianFitResult.mat'));

end
