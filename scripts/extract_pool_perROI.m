function pool_table = extract_pool_perROI(popt, ROI, pool_names, slice)
% extract_pool_perROI  Extract Lorentzian fit parameters for selected pools from an ROI
%
% USAGE:
%   pool_table = extract_pool_ROI(popt, ROI, pool_names)
%   pool_table = extract_pool_ROI(popt, ROI, pool_names, slice)
%
% INPUTS:
%   popt        - Fitted parameter array [X x Y x Z x nparams]
%   ROI         - Binary/logical mask [X x Y] or [X x Y x Z]
%   pool_names  - Cell array of pool names to extract, chosen from:
%                 'water', 'amide', 'NOE', 'MT', 'amine', 'glucose'
%   slice       - (optional) Slice index, default = 1
%
% OUTPUT:
%   pool_table  - MATLAB table with Amp, LW, Offset columns per requested pool
%
% EXAMPLES:
%   % Extract only glucose:
%   T = extract_pool_ROI(popt, ROI, {'glucose'});
%
%   % Extract glucose and MT:
%   T = extract_pool_ROI(popt, ROI, {'glucose','MT'});
%
%   % Extract all pools:
%   T = extract_pool_ROI(popt, ROI, {'water','amide','NOE','MT','amine','glucose'});

    %% --- Pool index map (pool k → popt indices 3k-1, 3k, 3k+1) ---
    pool_map = struct(...
        'water',   1, ...
        'amide',   2, ...
        'NOE',     3, ...
        'MT',      4, ...
        'amine',   5, ...
        'glucose', 6  ...
    );

    %% --- Defaults ---
    if nargin < 4
        slice = 1;
    end
    if nargin < 3 || isempty(pool_names)
        pool_names = {'glucose'};   % default: only glucose
    end
    if ischar(pool_names)
        pool_names = {pool_names};  % allow single string input
    end

    %% --- Build 2D ROI mask ---
    ROI_mask = logical(squeeze(ROI(:,:,min(slice,size(ROI,3)))));

    %% --- Extract requested pools ---
    col_names = {};
    data_matrix = [];

    for pp = 1:numel(pool_names)
        name = lower(pool_names{pp});

        if ~isfield(pool_map, name)
            warning('Pool "%s" not recognized. Skipping. Valid: water,amide,NOE,MT,amine,glucose', name);
            continue;
        end

        k     = pool_map.(name);
        i_amp = 3*k - 1;   % Amplitude index
        i_lw  = 3*k;       % Linewidth index
        i_dw  = 3*k + 1;   % Offset index

        amp_map = squeeze(popt(:,:,slice, i_amp));
        lw_map  = squeeze(popt(:,:,slice, i_lw));
        dw_map  = squeeze(popt(:,:,slice, i_dw));

        amp_ROI = amp_map(ROI_mask);
        lw_ROI  = lw_map(ROI_mask);
        dw_ROI  = dw_map(ROI_mask);

        data_matrix = [data_matrix, amp_ROI, lw_ROI, dw_ROI];

        col_names{end+1} = sprintf('%s_Amp',    name);
        col_names{end+1} = sprintf('%s_LW_ppm', name);
        col_names{end+1} = sprintf('%s_dw_ppm', name);
    end

    %% --- Remove voxels where ANY extracted value is NaN ---
    valid = all(~isnan(data_matrix), 2);
    data_matrix = data_matrix(valid, :);

    %% --- Build output table ---
    pool_table = array2table(data_matrix, 'VariableNames', col_names);

    %% --- Print summary ---
    fprintf('\n=== ROI Pool Summary: %d valid voxels ===\n', sum(valid));
    for cc = 1:numel(col_names)
        vals = data_matrix(:, cc);
        fprintf('  %-22s  mean = %8.4f   std = %8.4f   median = %8.4f\n', ...
            col_names{cc}, mean(vals), std(vals), median(vals));
    end

    %% --- Optional: save CSV ---
    fname = sprintf('pool_ROI_%s.csv', strjoin(pool_names, '_'));
    writetable(pool_table, fname);
    fprintf('\nSaved to %s\n', fname);

end
