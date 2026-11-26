function out = computeCESTRescueQC(dataDir, varargin)
% Rescue voxels excluded only by outer-offset problems + neighbourhood likelihood.
% Needs LorentzianFitResult.mat, QC_result.mat, QC_extra.mat. No refit. R2018a ok.
% Writes QC_final.mat, QC_final_summary.txt, ROI_<label>_QC2.mat (originals untouched).
def = CESTQCConfig();
p = inputParser;
addParameter(p,'RoiFile',def.RoiFile);
addParameter(p,'NTrim',2);
addParameter(p,'GlucWin',[0.5 3.0]);   % |ppm| window; CHECK against your glucose offsets
addParameter(p,'GlucResidMax',0.05);
addParameter(p,'RatioThresh',10);
addParameter(p,'GlucFloor',0.005);     % same units as gluc_integral_fullMap
addParameter(p,'SupportFactor',3);
addParameter(p,'MinSupport',2);        % other similar voxels needed in 5x5
addParameter(p,'MinNeighbors',3);
addParameter(p,'Suffix','_QC2');
parse(p,varargin{:}); o = p.Results;
roiFiles = cellstr(o.RoiFile);

F = load(fullfile(dataDir,'LorentzianFitResult.mat'),'popt','P','Segment','Z_corrExt','gluc_integral_fullMap');
Q = load(fullfile(dataDir,'QC_result.mat'),'badPixelMask','qcFlags');
E = load(fullfile(dataDir,'QC_extra.mat'),'extra'); ex = E.extra;
W = double(F.P.SEQ.w(:)); Z = double(F.Z_corrExt); [nx,ny,~,nW] = size(Z);
if numel(W) ~= nW, error('computeCESTRescueQC:offsetMismatch','%d vs %d in %s',numel(W),nW,dataDir); end
th = Q.qcFlags.thresholds;
tissue = (F.Segment(:,:,1) == 1);
basicBad = logical(Q.badPixelMask); finalBad = logical(ex.finalBad);
inner = (o.NTrim+1):(nW-o.NTrim);
gw = abs(W) >= o.GlucWin(1) & abs(W) <= o.GlucWin(2);

%% ---- A) rescue ----
cand = find(basicBad & ~logical(ex.extraBad) & tissue & ...
    ~Q.qcFlags.masks.Bound & ~Q.qcFlags.masks.Spatial);
rescued = false(nx,ny); rmseIn = nan(nx,ny); r2In = nan(nx,ny); gMax = nan(nx,ny);
for c = 1:numel(cand)
    [xi,yi] = ind2sub([nx ny], cand(c));
    z = reshape(Z(xi,yi,1,:),[],1); pp = double(squeeze(F.popt(xi,yi,1,:))); pp = pp(:);
    if any(~isfinite(z)) || any(~isfinite(pp)), continue; end
    [f,~] = fitmodelfunc_ANA(pp, F.P);
    yh = reshape(f(W),[],1); if numel(yh) ~= nW, continue; end
    r = z - yh; ri = r(inner); zi = z(inner);
    sst = sum((zi-mean(zi)).^2); if sst <= 0, continue; end
    rmseIn(xi,yi) = sqrt(mean(ri.^2));
    r2In(xi,yi)   = 1 - sum(ri.^2)/sst;
    gMax(xi,yi)   = max(abs(r(gw)));
    rescued(xi,yi) = min(zi) <= th.ZminThresh && max(zi) <= th.ZmaxThresh && ...
        rmseIn(xi,yi) <= th.RMSEThreshUsed && r2In(xi,yi) >= th.R2Thresh && ...
        gMax(xi,yi) <= o.GlucResidMax;
end

%% ---- ROIs ----
ROIs = cell(1,numel(roiFiles)); labels = cell(1,numel(roiFiles)); anyROI = false(nx,ny);
for k = 1:numel(roiFiles)
    [~,lab] = fileparts(roiFiles{k}); labels{k} = strrep(lab,'ROI_','');
    rp = fullfile(dataDir, roiFiles{k});
    if ~isfile(rp), continue; end
    S = load(rp);
    if isfield(S,'ROI') && isequal(size(S.ROI),[nx ny]), ROIs{k} = logical(S.ROI); anyROI = anyROI | ROIs{k}; end
end

%% ---- B) neighbourhood likelihood ----
G = double(F.gluc_integral_fullMap); if ndims(G) > 2, G = G(:,:,1); end
pass = ((tissue & ~finalBad) | rescued) & isfinite(G);
Gp = nan(nx+4,ny+4); Gp(3:end-2,3:end-2) = G;
Pp = false(nx+4,ny+4); Pp(3:end-2,3:end-2) = pass;
spatialBad = false(nx,ny); hotKeep = false(nx,ny); ratioMap = nan(nx,ny);
for idx = find(pass & anyROI)'
    [i,j] = ind2sub([nx ny], idx);
    blk = Gp(i:i+4, j:j+4); pb = Pp(i:i+4, j:j+4); pb(3,3) = false;
    n3 = blk(2:4,2:4); p3 = pb(2:4,2:4); v3 = n3(p3);
    if numel(v3) < o.MinNeighbors, continue; end
    v = max(G(i,j), o.GlucFloor); m = max(median(v3), o.GlucFloor);
    ratioMap(i,j) = v/m;
    if v/m >= o.RatioThresh || m/v >= o.RatioThresh
        rr = max(blk(pb), o.GlucFloor) / v;
        nSup = sum(rr <= o.SupportFactor & rr >= 1/o.SupportFactor);
        if nSup >= o.MinSupport, hotKeep(i,j) = true; else, spatialBad(i,j) = true; end
    end
end
finalBad2 = (finalBad & ~rescued) | spatialBad;

%% ---- per-ROI summary + ROI files ----
lines = {sprintf('CEST rescue+spatial QC | %s | nOffsets=%d | %s', dataDir, nW, datestr(now))};
lines{end+1} = sprintf('Rescue: NTrim=%d GlucWin=[%g %g] GlucResidMax=%g | Spatial: ratio>=%g, floor=%g, support>=%d within x%g (5x5), minNb=%d', ...
    o.NTrim,o.GlucWin(1),o.GlucWin(2),o.GlucResidMax,o.RatioThresh,o.GlucFloor,o.MinSupport,o.SupportFactor,o.MinNeighbors);
roi = struct('label',labels,'nROI',NaN,'nOldBad',NaN,'nRescued',NaN,'nSpatial',NaN,'nHotKept',NaN,'nClean',NaN,'mean',NaN,'sd',NaN);
newROI = cell(1,numel(roiFiles));
for k = 1:numel(roiFiles)
    if isempty(ROIs{k}), continue; end
    R = ROIs{k}; newROI{k} = R & ~finalBad2; gv = G(newROI{k});
    roi(k).nROI = nnz(R); roi(k).nOldBad = nnz(R & finalBad); roi(k).nRescued = nnz(R & rescued);
    roi(k).nSpatial = nnz(R & spatialBad); roi(k).nHotKept = nnz(R & hotKeep);
    roi(k).nClean = nnz(newROI{k}); roi(k).mean = mean(gv,'omitnan'); roi(k).sd = std(gv,'omitnan');
    lines{end+1} = sprintf('ROI %s: n %d | oldBad %d | rescued %d | spatialBad %d | hotspotKept %d | nClean %d (%.1f%% kept) | mean %.4f SD %.4f', ...
        labels{k}, roi(k).nROI, roi(k).nOldBad, roi(k).nRescued, roi(k).nSpatial, roi(k).nHotKept, ...
        roi(k).nClean, 100*roi(k).nClean/max(roi(k).nROI,1), roi(k).mean, roi(k).sd);
    [~,base] = fileparts(roiFiles{k});
    ROI = cast(newROI{k}, class(R)); %#ok<NASGU>
    save(fullfile(dataDir,[base o.Suffix '.mat']),'ROI');
end
cellfun(@(s) fprintf('%s\n',s), lines);
fid = fopen(fullfile(dataDir,'QC_final_summary.txt'),'w'); fprintf(fid,'%s\n',lines{:}); fclose(fid);
out = struct('finalBad2',finalBad2,'rescued',rescued,'spatialBad',spatialBad,'hotKeep',hotKeep, ...
    'ratioMap',ratioMap,'rmseInner',rmseIn,'r2Inner',r2In,'glucWinMaxResid',gMax,'roi',roi,'opt',o);
save(fullfile(dataDir,'QC_final.mat'),'out','finalBad2','-v7');
end