function extra = computeCESTExtraQC(dataDir, varargin)
% Extra no-refit QC on top of computeCESTQualityMetrics (needs QC_result.mat).
% Reads LorentzianFitResult.mat + QC_result.mat + ROI files. R2018a compatible.
% Writes: QC_extra.mat, QC_extra_summary.txt, ROI_<label>_QC.mat (originals untouched).
%
% Design notes:
%  - Resid and Phys use INTERIOR offsets only (3:end-2), so the two outermost
%    points on each side (+-5 and +-4.69 ppm; likely dummy-scan/steady-state
%    affected) cannot trigger them. Their problems are reported separately as
%    endBad (information only, NOT used for exclusion).
%  - Asym (endpoint asymmetry vs ROI median) is information only.
%  - Exclusion = Invalid  OR  (Resid + Phys + Corr >= MinFlagsSoft)
%  - Invalid also covers ROI voxels outside the tissue mask / NaN data (never evaluated by basic QC).
%  - extra.roi(k): per-ROI counts + mean/SD of gluc_integral_fullMap (saved fit) over the final clean ROI.
p = inputParser;
addParameter(p, 'RoiFile', {'ROI_ctx.mat','ROI_medulla.mat'});
addParameter(p, 'ResidKSigma', 4);     % max|resid| / robust sigma (interior offsets)
addParameter(p, 'ResidAbsMin', 0.05);  % ...AND max|resid| must exceed this (Z units)
addParameter(p, 'SigmaFloor', 0.01);   % floor for robust sigma
addParameter(p, 'AsymMADk', 4);        % endpoint asymmetry robust z within ROI (info only)
addParameter(p, 'EndResidAbs', 0.10);  % outer-point |resid| above this -> endBad (info only)
addParameter(p, 'NTrim', 2);           % outer points trimmed on EACH side for Resid/Phys
addParameter(p, 'CorrMin', 0.95);      % corr with ROI median spectrum
addParameter(p, 'ZPhysLo', -0.02);
addParameter(p, 'ZPhysHi', 1.10);
addParameter(p, 'MinFlagsSoft', 2);    % soft flags needed to exclude
addParameter(p, 'SaveROI', true);
addParameter(p, 'ROISuffix', '_QC');
parse(p, varargin{:}); o = p.Results;
roiFiles = cellstr(o.RoiFile);

F = load(fullfile(dataDir,'LorentzianFitResult.mat'), 'popt','P','Segment','Z_corrExt','gluc_integral_fullMap');
Q = load(fullfile(dataDir,'QC_result.mat'), 'badPixelMask');
W = double(F.P.SEQ.w(:)); Z = F.Z_corrExt;
[nx,ny,~,nW] = size(Z);
if numel(W) ~= nW
    error('computeCESTExtraQC:offsetMismatch','%d offsets in P.SEQ.w vs %d in Z_corrExt (%s).', numel(W), nW, dataDir);
end
if any(diff(W) <= 0), error('computeCESTExtraQC:unsorted','W not ascending in %s', dataDir); end
nt = o.NTrim;
if nW < 2*nt + 5, error('computeCESTExtraQC:tooFewOffsets','Only %d offsets in %s', nW, dataDir); end
tissue = (F.Segment(:,:,1) == 1);
inner  = (nt+1):(nW-nt);
outer  = [1:nt, (nW-nt+1):nW];

% ---- load ROIs ----
ROIs = cell(1,numel(roiFiles)); labels = cell(1,numel(roiFiles)); anyROI = false(nx,ny);
for k = 1:numel(roiFiles)
    [~,lab] = fileparts(roiFiles{k}); labels{k} = strrep(lab,'ROI_','');
    rp = fullfile(dataDir, roiFiles{k});
    if ~isfile(rp), warning('ROI not found: %s', rp); continue; end
    S = load(rp);
    if ~isfield(S,'ROI') || ~isequal(size(S.ROI), [nx ny])
        warning('Bad ROI variable/size in %s', rp); continue;
    end
    ROIs{k} = S.ROI; anyROI = anyROI | logical(S.ROI);
end

% ---- per-voxel metrics ----
valid = false(nx,ny); maxResid = nan(nx,ny); maxAbs = nan(nx,ny);
asym = nan(nx,ny); physBad = false(nx,ny); endBad = false(nx,ny);
resid3 = nan(nx,ny,nW); sigmaMap = nan(nx,ny);
for i = find(tissue)'
    [xi,yi] = ind2sub([nx ny], i);
    z = double(reshape(Z(xi,yi,1,:), [], 1));
    pp = double(squeeze(F.popt(xi,yi,1,:))); pp = pp(:);
    if any(~isfinite(z)) || any(~isfinite(pp)), continue; end
    [f,~] = fitmodelfunc_ANA(pp, F.P);
    yh = reshape(f(W), [], 1);
    if numel(yh) ~= nW, continue; end
    r = z - yh; ri = r(inner);
    s = max(1.4826*median(abs(ri - median(ri))), o.SigmaFloor);
    valid(xi,yi) = true;
    resid3(xi,yi,:) = r; sigmaMap(xi,yi) = s;
    maxAbs(xi,yi) = max(abs(ri)); maxResid(xi,yi) = maxAbs(xi,yi)/s;
    asym(xi,yi) = z(end) - z(1);
    zi = z(inner);
    physBad(xi,yi) = any(zi < o.ZPhysLo) || any(zi > o.ZPhysHi);
    endBad(xi,yi) = max(abs(r(outer))) > o.EndResidAbs;
end

flag.Invalid = (tissue | anyROI) & ~valid;
flag.Resid   = valid & (maxResid > o.ResidKSigma) & (maxAbs > o.ResidAbsMin);
flag.Phys    = valid & physBad;
flag.Asym    = false(nx,ny);   % info only
flag.Corr    = false(nx,ny);
flag.EndBad  = valid & endBad; % info only
corrMap = nan(nx,ny);
offBadFrac = nan(nW, numel(roiFiles));

for k = 1:numel(roiFiles)
    if isempty(ROIs{k}), continue; end
    R = logical(ROIs{k}) & valid; idx = find(R);
    if numel(idx) < 5, continue; end
    a = asym(R); m = median(a); s = max(1.4826*median(abs(a-m)), 0.01);
    flag.Asym(R) = flag.Asym(R) | (abs(a - m) > o.AsymMADk * s);
    A = nan(numel(idx), nW);
    for j = 1:numel(idx)
        [xi,yi] = ind2sub([nx ny], idx(j));
        A(j,:) = reshape(Z(xi,yi,1,:), 1, []);
    end
    med = median(A(:,inner),1);
    for j = 1:numel(idx)
        c = corrcoef(A(j,inner), med);
        corrMap(idx(j)) = c(1,2);
    end
    flag.Corr(idx) = flag.Corr(idx) | (corrMap(idx) < o.CorrMin) | isnan(corrMap(idx));
    for w = 1:nW
        rw = resid3(:,:,w);
        offBadFrac(w,k) = mean(abs(rw(R)) > o.ResidKSigma * sigmaMap(R) & abs(rw(R)) > o.ResidAbsMin);
    end
end

soft = double(flag.Resid) + double(flag.Phys) + double(flag.Corr);
extraBad = flag.Invalid | (soft >= o.MinFlagsSoft);
finalBad = logical(Q.badPixelMask) | extraBad;

% ---- per-ROI summary + new ROI files ----
lines = {sprintf('CEST extra QC | %s | N offsets = %d | %s', dataDir, nW, datestr(now))};
lines{end+1} = sprintf(['Params: ResidK=%g ResidAbsMin=%g CorrMin=%g softFlags>=%d NTrim=%d | ' ...
    'Resid/Phys on interior offsets | Asym and EndBad are info only'], ...
    o.ResidKSigma, o.ResidAbsMin, o.CorrMin, o.MinFlagsSoft, nt);
newROI = cell(1,numel(roiFiles));
roiStats = struct('label',labels,'nROI',NaN,'nInvalid',NaN,'nResid',NaN,'nPhys',NaN,'nCorr',NaN, ...
    'nExtraOnly',NaN,'nFinalBad',NaN,'nClean',NaN,'meanFinal',NaN,'sdFinal',NaN);
for k = 1:numel(roiFiles)
    if isempty(ROIs{k}), continue; end
    R = logical(ROIs{k});
    nROI = nnz(R);
    newROI{k} = R & ~finalBad;
    gv = F.gluc_integral_fullMap(newROI{k});
    roiStats(k).nROI = nROI;
    roiStats(k).nInvalid = nnz(R&flag.Invalid);
    roiStats(k).nResid = nnz(R&flag.Resid);
    roiStats(k).nPhys = nnz(R&flag.Phys);
    roiStats(k).nCorr = nnz(R&flag.Corr);
    roiStats(k).nExtraOnly = nnz(R & extraBad & ~logical(Q.badPixelMask));
    roiStats(k).nFinalBad = nnz(R&finalBad);
    roiStats(k).nClean = nnz(newROI{k});
    roiStats(k).meanFinal = mean(gv(:),'omitnan');
    roiStats(k).sdFinal = std(gv(:),'omitnan');
    lines{end+1} = sprintf(['ROI %s: nROI %d | Invalid %d | Resid %d | Phys %d | Corr %d | ' ...
        'Asym(info) %d | EndBad(info) %d | oldBad %d | extraBad %d | finalBad %d | nClean %d (%.1f%% kept)'], ...
        labels{k}, nROI, nnz(R&flag.Invalid), nnz(R&flag.Resid), nnz(R&flag.Phys), nnz(R&flag.Corr), ...
        nnz(R&flag.Asym), nnz(R&flag.EndBad), nnz(R&Q.badPixelMask), nnz(R&extraBad), ...
        nnz(R&finalBad), nnz(newROI{k}), 100*nnz(newROI{k})/max(nROI,1));
    if o.SaveROI
        [~,base] = fileparts(roiFiles{k});
        ROI = cast(newROI{k}, class(ROIs{k})); %#ok<NASGU>
        save(fullfile(dataDir, [base o.ROISuffix '.mat']), 'ROI');
    end
end
cellfun(@(s) fprintf('%s\n', s), lines);
fid = fopen(fullfile(dataDir,'QC_extra_summary.txt'),'w');
fprintf(fid,'%s\n', lines{:}); fclose(fid);

extra = struct('flag',flag,'maxResidSigma',maxResid,'maxAbsResid',maxAbs,'asymMap',asym, ...
    'corrMap',corrMap,'endBad',endBad,'offBadFrac',offBadFrac,'W',W,'extraBad',extraBad, ...
    'finalBad',finalBad,'newROI',{newROI},'roiLabels',{labels},'roi',roiStats,'opt',o);
save(fullfile(dataDir,'QC_extra.mat'), 'extra', '-v7');
end