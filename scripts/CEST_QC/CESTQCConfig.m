function qc = CESTQCConfig()
% CESTQCCONFIG  Single source of truth for ALL CEST QC thresholds.
%
%   qc = CESTQCConfig()
%
%   EDIT THRESHOLDS HERE AND ONLY HERE.
%   computeCESTQualityMetrics reads its defaults from this file, so
%   WZ_QueueCEST_QC.m (batch) and plotCESTQCBeforeAfter (interactive)
%   automatically use the same numbers. Name-value overrides passed to
%   computeCESTQualityMetrics / QCOptions still win, but they are
%   recorded as "overrides" in QC_summary.txt so nothing is hidden.
%
%   Field names below are IDENTICAL to the name-value parameter names of
%   computeCESTQualityMetrics. Do not rename.
%
%   Criterion directions (as implemented in computeCESTQualityMetrics):
%     Zmin    : flag if  Zmin_map  >  ZminThresh
%     Zmax    : flag if  Zmax_map  >  ZmaxThresh
%     RMSE    : flag if  RMSE_map  >  RMSEThreshUsed
%               (RMSEThreshUsed = prctile(RMSE_map(tissueMask), RMSEPercentile)
%                if UseRMSEPercentile, otherwise = RMSEThresh)
%     R2      : flag if  R2_map    <  R2Thresh
%     Bound   : flag if  any(|popt - lb| < BoundTol) or any(|popt - ub| < BoundTol)
%     Spatial : flag if  |g - medfilt2(g)| > SpatialMADMultiplier * medfilt2(|g - medfilt2(g)|)
%               with g = gluc_integral_fullMap and window = SpatialWindow
%
%   MATLAB R2018a compatible.

% ---------------- ROI files evaluated per session ----------------
qc.RoiFile              = {'ROI_ctx.mat', 'ROI_medulla.mat'};

% ---------------- Z-spectrum range criteria ----------------
qc.ZminThresh           = 0.1;      % flag if Zmin_map > ZminThresh
qc.ZmaxThresh           = 1.2;      % flag if Zmax_map > ZmaxThresh

% ---------------- residual criteria ----------------
qc.UseRMSEPercentile    = false;     % true : data-relative (percentile) RMSE cut
                                    % false: absolute RMSE cut (RMSEThresh)
qc.RMSEPercentile       = 70;       % used only if UseRMSEPercentile == true
qc.RMSEThresh           = 0.05;      % used only if UseRMSEPercentile == false
qc.R2Thresh             = 0.95;      % flag if R2_map < R2Thresh

% ---------------- fit-parameter-at-bound criterion ----------------
qc.UseBoundCheck        = true;
qc.BoundTol             = 1e-14;

% ---------------- spatial outlier criterion ----------------
qc.UseSpatialOutlier    = false;
qc.SpatialMADMultiplier = 5;
qc.SpatialWindow        = [3 3];    % medfilt2 window (odd sizes only)

end
