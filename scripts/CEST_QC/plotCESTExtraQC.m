function plotCESTExtraQC(dataDir, varargin)
% Visual check of computeCESTExtraQC output. No refit. R2018a compatible.
p = inputParser;
addParameter(p, 'MaxSpectra', 40);     % max spectra drawn per group per ROI
addParameter(p, 'Save', true);
addParameter(p, 'RotN', 3);
addParameter(p, 'Visible', 'on');
parse(p, varargin{:}); o = p.Results;

E = load(fullfile(dataDir,'QC_extra.mat'),'extra'); extra = E.extra;
Q = load(fullfile(dataDir,'QC_result.mat'),'badPixelMask');
F = load(fullfile(dataDir,'LorentzianFitResult.mat'),'M0_stack','Segment','Z_corrExt');
if isfield(F,'M0_stack'), M0 = double(F.M0_stack); else, M0 = []; end
Z = F.Z_corrExt; [nx,ny,~,nW] = size(Z);
W = extra.W; opt = extra.opt;
oldBad = logical(Q.badPixelMask); newBad = extra.extraBad & ~oldBad; fin = extra.finalBad;
labels = extra.roiLabels; nR = numel(labels);
[~,tag] = fileparts(dataDir); if isempty(tag), [~,tag] = fileparts(fileparts(dataDir)); end

% load ROIs (original ones)
ROIs = cell(1,nR);
for k = 1:nR
    rp = fullfile(dataDir, ['ROI_' labels{k} '.mat']);
    if isfile(rp), S = load(rp); ROIs{k} = logical(S.ROI); end
end
bg = M0; if isempty(bg) || ~isequal(size(bg),[nx ny]), bg = mean(squeeze(Z(:,:,1,:)),3); end
clim = [min(bg(:)) prctile(bg(:),99)]; if clim(2) <= clim(1), clim = [0 1]; end
cols = lines(max(nR,3));

% ---------- rotated copies for map display (spectra stay in raw frame) ----------
bgR = rot90(bg,o.RotN); oldBadR = rot90(oldBad,o.RotN); newBadR = rot90(newBad,o.RotN);
[nxr,nyr] = size(bgR);
ROIsR = cell(1,nR);
for k = 1:nR, if ~isempty(ROIs{k}), ROIsR{k} = rot90(ROIs{k},o.RotN); end, end
solid = @(rgb) cat(3, rgb(1)*ones(nxr,nyr), rgb(2)*ones(nxr,nyr), rgb(3)*ones(nxr,nyr));

% ---------- Fig 1: exclusion map ----------
f1 = figure('Visible',o.Visible,'Position',[50 50 600 560]);
imagesc(bgR, clim); colormap(gca, gray); axis image off; hold on;
h1 = image(gca, solid([1 0 0]));   set(h1,'AlphaData',0.55*double(oldBadR));
h2 = image(gca, solid([1 0.6 0])); set(h2,'AlphaData',0.55*double(newBadR));
for k = 1:nR
    if isempty(ROIsR{k}), continue; end
    B = bwboundaries(ROIsR{k});
    for b = 1:numel(B), plot(B{b}(:,2),B{b}(:,1),'-','Color',cols(k,:),'LineWidth',1.5); end
end
title(sprintf('%s\nred=old QC (%d) | orange=extra QC only (%d) | ROI outlines: %s', tag, ...
    nnz(oldBad), nnz(newBad), strjoin(labels,'/')), 'Interpreter','none','FontSize',10);

% ---------- Fig 2: per-flag maps ----------
names = {'Invalid','Resid','Asym','Phys','Corr'};
f2 = figure('Visible',o.Visible,'Position',[50 50 1200 420]);
for i = 1:numel(names)
    subplot(1,5,i); imagesc(bgR, clim); colormap(gca, gray); axis image off; hold on;
    m = extra.flag.(names{i});
    h = image(gca, solid([1 0 0])); set(h,'AlphaData',0.7*double(rot90(m,o.RotN)));
    for k = 1:nR
        if isempty(ROIsR{k}), continue; end
        B = bwboundaries(ROIsR{k});
        for b = 1:numel(B), plot(B{b}(:,2),B{b}(:,1),'-','Color',cols(k,:),'LineWidth',1); end
    end
    title(sprintf('%s (%d)', names{i}, nnz(m)),'FontSize',10);
end

% ---------- Fig 3: spectra per ROI ----------
f3 = figure('Visible',o.Visible,'Position',[50 50 450*nR 450]);
for k = 1:nR
    subplot(1,nR,k); hold on; box on;
    if isempty(ROIs{k}), title([labels{k} ' (missing)']); continue; end
    R = ROIs{k};
    keepIdx = find(R & ~fin); badIdx = find(R & fin);
    keepIdx = keepIdx(randperm(numel(keepIdx), min(numel(keepIdx), o.MaxSpectra)));
    badIdx  = badIdx(randperm(numel(badIdx),  min(numel(badIdx),  o.MaxSpectra)));
    for j = 1:numel(keepIdx)
        [xi,yi] = ind2sub([nx ny], keepIdx(j));
        plot(W, squeeze(Z(xi,yi,1,:)), '-', 'Color', [0.6 0.6 0.6]);
    end
    for j = 1:numel(badIdx)
        [xi,yi] = ind2sub([nx ny], badIdx(j));
        c = [1 0 0]; if newBad(xi,yi), c = [1 0.6 0]; end
        plot(W, squeeze(Z(xi,yi,1,:)), '-', 'Color', c);
    end
    kk = R & ~fin;
    if any(kk(:))
        A = nan(nnz(kk), nW); ii = find(kk);
        for j = 1:numel(ii), [xi,yi] = ind2sub([nx ny], ii(j)); A(j,:) = squeeze(Z(xi,yi,1,:))'; end
        plot(W, median(A,1,'omitnan'), 'k-', 'LineWidth', 2.5);
    end
    set(gca,'XDir','reverse'); xlabel('offset (ppm)'); ylabel('Z');
    title(sprintf('%s: kept %d, excluded %d (grey=kept, red=old, orange=extra, black=median kept)', ...
        labels{k}, nnz(R&~fin), nnz(R&fin)), 'Interpreter','none','FontSize',8);
end

% ---------- Fig 4: offset diagnosis ----------
f4 = figure('Visible',o.Visible,'Position',[50 50 600 380]); hold on; box on;
for k = 1:nR, plot(W, extra.offBadFrac(:,k), 'o-', 'Color', cols(k,:), 'LineWidth', 1.5); end
set(gca,'XDir','reverse'); xlabel('offset (ppm)'); ylabel('fraction of ROI voxels with outlier residual');
legend(labels,'Location','best'); title('Offset-wise outlier fraction (frame-level motion check)');

% ---------- Fig 5: metric distributions ----------
f5 = figure('Visible',o.Visible,'Position',[50 50 1100 350]);
tis = ~isnan(extra.maxResidSigma);
subplot(1,3,1); histogram(extra.maxResidSigma(tis), 30); hold on; yl = ylim;
plot(opt.ResidKSigma*[1 1], yl, 'r--'); xlabel('max|resid| / sigma'); ylabel('voxels');
subplot(1,3,2); histogram(extra.asymMap(tis), 30); xlabel('Z(+end) - Z(-end)'); ylabel('voxels');
subplot(1,3,3); histogram(extra.corrMap(~isnan(extra.corrMap)), 30); hold on; yl = ylim;
plot(opt.CorrMin*[1 1], yl, 'r--'); xlabel('corr with ROI median spectrum'); ylabel('voxels');

if o.Save
    figs = {f1,f2,f3,f4,f5}; fn = {'exclusion_map','flag_maps','spectra','offset_diag','metric_hist'};
    for i = 1:numel(figs)
        set(figs{i},'PaperPositionMode','auto');
        print(figs{i}, fullfile(dataDir, ['QC_extra_' fn{i} '.png']), '-dpng', '-r150');
    end
end
if strcmpi(o.Visible,'off'), close([f1 f2 f3 f4 f5]); end
end
