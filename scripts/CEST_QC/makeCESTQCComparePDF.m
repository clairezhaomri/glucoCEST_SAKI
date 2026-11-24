function makeCESTQCComparePDF()
% Before QC | old QC | extra QC v1 | extra QC v2, overlaid on M0 image, one page per session.
% R2018a compatible. No refit. Writes a new timestamped folder.
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
Filter = 'Ctrl'; Scale = 4; Margin = 5;

d = dir(fullfile(Basepath, '**', 'QC_result.mat'));
dirs = sort(unique({d.folder})); dirs = dirs(contains(dirs, Filter));
% dirs = dirs(1:2);   % <- uncomment to test on 2 sessions first
fprintf('%d sessions\n', numel(dirs));
ReportBase = fullfile(Basepath, 'QC_Report'); if ~isfolder(ReportBase), mkdir(ReportBase); end
RunDir = fullfile(ReportBase, ['QCcompare_' datestr(now,'yyyymmdd_HHMMSS')]); mkdir(RunDir);
pageDir = fullfile(RunDir, 'pages'); mkdir(pageDir);
PDF = fullfile(RunDir, 'QC_compare_before_old_v1_v2.pdf');

rows = cell(0,10); pageFiles = {}; page = 0;
for i = 1:numel(dirs)
    f = [];
    try
        S = computeSession(dirs{i});
        f = renderPage(S, strrep(dirs{i}, Basepath, ''), Scale, Margin);
        pf = fullfile(pageDir, sprintf('page_%03d.pdf', page+1));
        savePagePDF(f, pf);
        close(f); page = page + 1; pageFiles{end+1} = pf; %#ok<AGROW>
        for k = 1:2
            gm = @(m) mean(S.g(m), 'omitnan');
            rows(end+1,:) = {strrep(dirs{i},Basepath,''), S.names{k}, page, ...
                nnz(S.M{1}{k}), nnz(S.M{2}{k}), nnz(S.M{3}{k}), nnz(S.M{4}{k}), ...
                gm(S.M{1}{k}), gm(S.M{4}{k}), S.diffSaved(k)}; %#ok<AGROW>
        end
        fprintf('[%d/%d] page %d: %s\n', i, numel(dirs), page, dirs{i});
    catch ME
        fprintf('FAILED %s: %s\n', dirs{i}, ME.message);
        if ~isempty(f) && isgraphics(f), close(f); end
    end
end
T = cell2table(rows, 'VariableNames', {'Session','ROI','Page','nBefore','nOldQC','nV1','nV2', ...
    'MeanBefore','MeanV2','V2DiffVsSavedROI'});
writetable(T, fullfile(RunDir, 'QC_compare_summary.csv'));

if page > 0 && mergePDFs(pageFiles, PDF)
    fprintf('\nDone: %d pages\n%s\n', page, PDF);
else
    fprintf('\n%d single-page PDFs written, automatic merge failed or no pages.\nMerge manually: %s\n', page, pageDir);
end
end

function savePagePDF(f, pf)
pos = get(f, 'Position'); w = pos(3)/100; h = pos(4)/100;
set(f, 'PaperUnits','inches', 'PaperSize',[w h], 'PaperPosition',[0 0 w h]);
print(f, pf, '-dpdf', '-r150');
end

function ok = mergePDFs(files, out)
setenv('PATH', [getenv('PATH') ':/usr/local/bin:/opt/homebrew/bin']);
q = @(c) ['"' c '"'];
args = strjoin(cellfun(q, files, 'UniformOutput', false), ' ');
ok = false;
joinTool = '/System/Library/Automator/Combine PDF Pages.action/Contents/MacOS/join';
if isfile(joinTool)
    system([q(joinTool) ' -o ' q(out) ' ' args]);
    ok = isfile(out);
end
if ~ok && system('which pdfunite > /dev/null 2>&1') == 0
    system(['pdfunite ' args ' ' q(out)]);
    ok = isfile(out);
end
if ~ok && system('which gs > /dev/null 2>&1') == 0
    system(['gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=' q(out) ' ' args]);
    ok = isfile(out);
end
end

function S = computeSession(dataDir)
F = load(fullfile(dataDir,'LorentzianFitResult.mat'));
Q = load(fullfile(dataDir,'QC_result.mat'), 'badPixelMask');
W = double(F.P.SEQ.w(:)); Z = double(F.Z_corrExt);
[nx,ny,~,nW] = size(Z); Z3 = reshape(Z, nx, ny, nW);
g = double(F.gluc_integral_fullMap); bad0 = logical(Q.badPixelMask);
tissue = (F.Segment(:,:,1) == 1);
bg = []; if isfield(F,'M0_stack'), bg = double(squeeze(F.M0_stack)); end
if isempty(bg) || ~isequal(size(bg), [nx ny]), bg = mean(Z3,3); end

names = {'ctx','medulla'}; orig = cell(1,2); saved2 = cell(1,2);
for k = 1:2
    orig{k}   = loadROI(fullfile(dataDir, ['ROI_' names{k} '.mat']), [nx ny]);
    saved2{k} = loadROI(fullfile(dataDir, ['ROI_' names{k} '_QC.mat']), [nx ny]);
end
anyROI = orig{1} | orig{2};
if ~any(anyROI(:)), error('No ROI found'); end

R3 = nan(nx,ny,nW); valid = false(nx,ny);
for i = find(tissue)'
    [xi,yi] = ind2sub([nx ny], i);
    z = squeeze(Z3(xi,yi,:)); z = z(:);
    pp = double(squeeze(F.popt(xi,yi,1,:))); pp = pp(:);
    if any(~isfinite(z)) || any(~isfinite(pp)), continue; end
    [fh,~] = fitmodelfunc_ANA(pp, F.P);
    yh = reshape(fh(W), [], 1);
    if numel(yh) ~= nW, continue; end
    R3(xi,yi,:) = reshape(z - yh, 1, 1, nW); valid(xi,yi) = true;
end
xb1 = extraBad(Z3, R3, valid, tissue, anyROI, orig, 1, false);
xb2 = extraBad(Z3, R3, valid, tissue, anyROI, orig, 2, true);

M = cell(1,4); diffSaved = nan(1,2);
for k = 1:2
    M{1}{k} = orig{k};
    M{2}{k} = orig{k} & ~bad0;
    M{3}{k} = orig{k} & ~(bad0 | xb1);
    M{4}{k} = orig{k} & ~(bad0 | xb2);
    if any(saved2{k}(:)), diffSaved(k) = nnz(M{4}{k} ~= saved2{k}); end
end
S = struct('g',g,'bg',bg,'orig',{orig},'M',{M},'names',{names},'diffSaved',diffSaved,'anyROI',anyROI,'Z3',Z3,'W',W);
end

function f = renderPage(S, label, ~, Margin)
RotN = 3; climBg = [0 150]; climG = [0 0.2]; alphaG = 0.7;
cmJ = jet(256);
[nx,ny] = size(S.g);
gR = rot90(S.g,RotN); bgR = rot90(S.bg,RotN); anyR = rot90(S.anyROI,RotN);
[nr,nc] = size(gR);
[r,c] = find(anyR);
r1 = max(min(r)-Margin,1); r2 = min(max(r)+Margin,nr);
c1 = max(min(c)-Margin,1); c2 = min(max(c)+Margin,nc);
W = S.W(:); nW = numel(W);

f = figure('Visible','off','Color','w','Units','pixels','Position',[50 50 2400 1800]);
annotation(f,'textbox',[.01 .955 .98 .04],'String',label,'Interpreter','none', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',14);
annotation(f,'textbox',[.01 .925 .98 .03],'EdgeColor','none','HorizontalAlignment','center','FontSize',10, ...
    'String',['Maps: glucose integral (a.u.), fixed 0-0.2, only voxels in each stage''s ROI colored. ' ...
    'Spectra: gray = voxel Z-spectra (max 80 shown), black = mean, band = +/-SD, colored = same 3 example voxels (10/50/90th pct glucose of final ROI).']);
sets = {'Before QC','After old QC','After extra QC v1','After extra QC v2'};
mapBottom = [0.69 0.215]; mapH = 0.21; spBottom = [0.50 0.055]; spH = 0.13;
exCols = [0.85 0.10 0.10; 0.10 0.30 0.85; 0.10 0.60 0.20];

for rr = 1:2
    if ~any(S.orig{rr}(:))
        ax = axes('Parent',f,'Position',[0.01 mapBottom(rr) 0.2 mapH]); axis(ax,'off');
        text(ax,0,.5,[S.names{rr} ': no ROI']); continue;
    end
    idxF = find(S.M{4}{rr} & ~isnan(S.g)); ex = [];
    if numel(idxF) >= 3
        [~,o] = sort(S.g(idxF));
        ii = max(round([0.1 0.5 0.9]*numel(idxF)), 1);
        ex = idxF(o(ii));
    end
    for s = 1:4
        mask = S.M{s}{rr};
        % ---- map ----
        ax = axes('Parent',f,'Position',[0.01+(s-1)*0.238, mapBottom(rr), 0.22, mapH]);
        bgc = bgR(r1:r2,c1:c2);
        base = min(max((bgc - climBg(1))/(climBg(2)-climBg(1)),0),1); RGB = repmat(base,[1 1 3]);
        gc = gR(r1:r2,c1:c2);
        mR = rot90(mask,RotN); show = mR(r1:r2,c1:c2) & ~isnan(gc);
        gn = min(max((gc - climG(1))/(climG(2)-climG(1)),0),1);
        idx = round(gn*255) + 1; idx(isnan(gc)) = 1;
        C = reshape(cmJ(idx(:),:), [size(gc) 3]);
        for ch = 1:3
            A = RGB(:,:,ch); B = C(:,:,ch);
            A(show) = (1-alphaG)*A(show) + alphaG*B(show);
            RGB(:,:,ch) = A;
        end
        image(ax, RGB); axis(ax,'image'); axis(ax,'off');
        gk = S.g(mask); gk = gk(~isnan(gk));
        title(ax, sprintf('%s | %s\nMean %.4f, SD %.4f (n=%d)', S.names{rr}, sets{s}, ...
            mean(gk), std(gk), numel(gk)), 'FontSize', 11, 'FontWeight','normal');

        % ---- spectra ----
        ax2 = axes('Parent',f,'Position',[0.035+(s-1)*0.238, spBottom(rr), 0.195, spH]);
        hold(ax2,'on'); box(ax2,'on');
        ids = find(mask); A = zeros(numel(ids), nW);
        for j = 1:numel(ids)
            [xi,yi] = ind2sub([nx ny], ids(j)); A(j,:) = squeeze(S.Z3(xi,yi,:))';
        end
        A = A(all(isfinite(A),2),:);
        if ~isempty(A)
            s0 = rng; rng(1); sel = randperm(size(A,1), min(80,size(A,1))); rng(s0);
            plot(ax2, W, A(sel,:)', '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
            m = mean(A,1); sd = std(A,0,1);
            fill([W; flipud(W)], [(m+sd)'; flipud((m-sd)')], [1 0.6 0.6], 'FaceAlpha',0.35, 'EdgeColor','none', 'Parent', ax2);
            plot(ax2, W, m, 'k-', 'LineWidth', 2);
            for e = 1:numel(ex)
                [xi,yi] = ind2sub([nx ny], ex(e));
                plot(ax2, W, squeeze(S.Z3(xi,yi,:)), '-o', 'Color', exCols(e,:), 'MarkerSize', 3, 'LineWidth', 1);
            end
            title(ax2, sprintf('Z-spectra (n=%d) | mean SD over offsets %.3f', size(A,1), mean(sd)), ...
                'FontSize', 10, 'FontWeight','normal');
        end
        set(ax2,'XDir','reverse','FontSize',9); xlim(ax2,[min(W) max(W)]); ylim(ax2,[0 1.1]); grid(ax2,'on');
        if rr == 2, xlabel(ax2,'Offset (ppm)'); end
        if s == 1, ylabel(ax2,'Z'); end
    end
    axc = axes('Parent',f,'Position',[0.955 mapBottom(rr)+0.02 0.008 0.17]);
    image(axc, [0 1], [climG(1) climG(2)], reshape(cmJ,[256 1 3]));
    set(axc, 'YDir','normal', 'XTick',[], 'YAxisLocation','right', 'FontSize', 9);
end
drawnow;
end

function xb = extraBad(Z3, R3, valid, tissue, anyROI, orig, nt, corrInner)
[nx,ny,nW] = size(Z3);
inner = (nt+1):(nW-nt); cidx = inner; if ~corrInner, cidx = 1:nW; end
Ri = R3(:,:,inner);
s  = max(1.4826*median(abs(Ri - median(Ri,3)),3), 0.01);
mx = max(abs(Ri),[],3);
fR = valid & (mx./s > 4) & (mx > 0.05);
Zi = Z3(:,:,inner);
fP = valid & any(Zi < -0.02 | Zi > 1.10, 3);
fC = false(nx,ny);
for k = 1:numel(orig)
    R = orig{k} & valid; idx = find(R);
    if numel(idx) < 5, continue; end
    A = zeros(numel(idx), nW);
    for j = 1:numel(idx)
        [xi,yi] = ind2sub([nx ny], idx(j)); A(j,:) = squeeze(Z3(xi,yi,:))';
    end
    medv = median(A(:,cidx),1); cr = nan(numel(idx),1);
    for j = 1:numel(idx)
        cc = corrcoef(A(j,cidx), medv); cr(j) = cc(1,2);
    end
    fC(idx(cr < 0.95)) = true;
end
inv  = (tissue | anyROI) & ~valid;
soft = double(fR) + double(fP) + double(fC);
xb = inv | (soft >= 2);
end

function R = loadROI(fn, sz)
R = false(sz);
if ~isfile(fn), return; end
S = load(fn);
if isfield(S,'ROI') && isequal(size(S.ROI), sz), R = logical(S.ROI); end
end