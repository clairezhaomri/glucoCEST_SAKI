function reviewCESTExtraQC(dataDir, maxPanels)
% Minimal visual review of voxels excluded ONLY by the extra QC. R2018a compatible.
if nargin < 2, maxPanels = 16; end
E = load(fullfile(dataDir,'QC_extra.mat'),'extra'); ex = E.extra;
Q = load(fullfile(dataDir,'QC_result.mat'),'badPixelMask');
F = load(fullfile(dataDir,'LorentzianFitResult.mat'),'popt','P','Z_corrExt');
W = ex.W(:); sz = size(Q.badPixelMask);

inROI = false(sz);
for k = 1:numel(ex.roiLabels)
    rp = fullfile(dataDir, ['ROI_' ex.roiLabels{k} '.mat']);
    if isfile(rp), S = load(rp); inROI = inROI | logical(S.ROI); end
end

nHi = nnz(inROI & (ex.asymMap > 5 | ex.asymMap < -5));
fprintf('ROI voxels with |Z(+end)-Z(-end)| > 5 (broken Z): %d\n', nHi);

idx = find(ex.extraBad & ~logical(Q.badPixelMask) & inROI);
fprintf('Voxels excluded ONLY by extra QC (in ROI): %d\n', numel(idx));
if isempty(idx), return; end
idx = idx(1:min(numel(idx), maxPanels));
n = numel(idx); nc = 4; nr = ceil(n/nc);

figure('Position',[50 50 1200 230*nr]);
for j = 1:n
    [xi,yi] = ind2sub(sz, idx(j));
    z  = double(squeeze(F.Z_corrExt(xi,yi,1,:))); z = z(:);
    pp = double(squeeze(F.popt(xi,yi,1,:)));      pp = pp(:);
    [f,~] = fitmodelfunc_ANA(pp, F.P); yh = f(W); yh = yh(:);
    subplot(nr,nc,j);
    plot(W,z,'ko',W,yh,'r-',W,z-yh-0.2,'b-'); set(gca,'XDir','reverse');
    names = {'Invalid','Resid','Asym','Phys','Corr'}; fl = '';
    for q = 1:numel(names)
        if ex.flag.(names{q})(xi,yi), fl = [fl names{q} ' ']; end
    end
    title(sprintf('(%d,%d) %s',xi,yi,fl),'FontSize',8);
    ylim([-0.5 1.3]);
end
set(gcf,'PaperPositionMode','auto');
print(gcf, fullfile(dataDir,'QC_extra_review.png'), '-dpng', '-r150');
end