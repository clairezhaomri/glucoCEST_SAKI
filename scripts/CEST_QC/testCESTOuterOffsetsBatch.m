function [All, Summary, outDir] = testCESTOuterOffsetsBatch(rootDir, varargin)
% MATLAB R2018a. Existing-fit diagnostic; NO refit or production-mask writes.
% rootDir: a CEST session, a subject folder, or the cohort Basepath.
% Requires CESTQCConfig, computeCESTQualityMetrics, fitmodelfunc_ANA.
p=inputParser;
addParameter(p,'RoiFile',{});
addParameter(p,'OutputRoot','');
parse(p,varargin{:}); o=p.Results;
cfg=CESTQCConfig();
if isempty(o.RoiFile), o.RoiFile=cfg.RoiFile; end
if ischar(o.RoiFile)||isstring(o.RoiFile), o.RoiFile=cellstr(o.RoiFile); end
if cfg.UseRMSEPercentile
    error('Absolute RMSE QC required. Percentile thresholds are not supported.');
end
if exist('fitmodelfunc_ANA','file')==0, error('fitmodelfunc_ANA missing from MATLAB path.'); end
sessions=local_sessions(rootDir);
if isempty(sessions), error('No CEST_base/CEST_infu session folders found.'); end
if isempty(o.OutputRoot), o.OutputRoot=rootDir; end
if ~exist(o.OutputRoot,'dir'), mkdir(o.OutputRoot); end
[~,token]=fileparts(tempname(o.OutputRoot));
outDir=fullfile(o.OutputRoot,['OuterOffsetTest_' datestr(now,'yyyymmdd_HHMMSS') '_' token]);
mkdir(outDir);
All=table(); Summary=table();
Log=cell(numel(sessions),3);
for s=1:numel(sessions)
    d=sessions{s}; Log{s,1}=d;
    fprintf('[%d/%d] %s\n',s,numel(sessions),d);
    try
        [T,U]=local_session(d,o.RoiFile);
        sessionOut=fullfile(outDir,sprintf('session_%04d',s)); mkdir(sessionOut);
        writetable(T,fullfile(sessionOut,'voxels.csv'));
        writetable(U,fullfile(sessionOut,'roi_summary.csv'));
        save(fullfile(sessionOut,'diagnostic.mat'),'T','U','d','cfg','-v7');
        if isempty(All), All=T; else, All=[All;T]; end %#ok<AGROW>
        if isempty(Summary), Summary=U; else, Summary=[Summary;U]; end %#ok<AGROW>
        Log{s,2}='OK'; Log{s,3}='';
    catch ME
        Log{s,2}='ERROR'; Log{s,3}=ME.message;
        warning('OuterOffsetTest:sessionFailed','%s: %s',d,ME.message);
    end
    writetable(cell2table(Log(1:s,:),'VariableNames',{'SessionPath','Status','Message'}),...
        fullfile(outDir,'session_log.csv'));
end
if ~isempty(All), writetable(All,fullfile(outDir,'ALL_voxels.csv')); end
if ~isempty(Summary), writetable(Summary,fullfile(outDir,'ALL_roi_summary.csv')); end
save(fullfile(outDir,'batch_settings.mat'),'rootDir','sessions','cfg','o','-v7');
fprintf('\nFinished: %d OK, %d failed.\n',sum(strcmp(Log(:,2),'OK')),sum(strcmp(Log(:,2),'ERROR')));
fprintf('Diagnostic only. Parameters, g maps and production QC masks unchanged.\n%s\n',outDir);
end

function [T,U]=local_session(d,roiFiles)
F=load(fullfile(d,'LorentzianFitResult.mat'),...
    'popt','P','Z_corrExt','Segment','gluc_integral_fullMap');
W=double(F.P.SEQ.w(:)); N=numel(W);
[nx,ny,~,np]=size(F.popt);
if N<4 || size(F.Z_corrExt,4)~=N, error('Invalid offset count.'); end
if size(F.popt,3)~=1 || size(F.Z_corrExt,3)~=1 || np<18
    error('Expected single-slice model with glucose parameters 17/18.');
end
if any(~isfinite(W))||numel(unique(W))~=N, error('Invalid/duplicate offsets.'); end
[negPPM,kn]=min(W); [posPPM,kp]=max(W);
if kn==kp, error('Endpoints coincide.'); end
masks=cell(numel(roiFiles),1); unionMask=false(nx,ny);
for k=1:numel(roiFiles)
    R=load(fullfile(d,roiFiles{k}),'ROI');
    if ~isfield(R,'ROI')||~isequal(size(R.ROI),[nx ny]), error('Invalid ROI: %s',roiFiles{k}); end
    masks{k}=isfinite(R.ROI)&R.ROI~=0;
    unionMask=unionMask|masks{k};
end
tissue=F.Segment(:,:,1)==1;
ids=find(unionMask & tissue);
if isempty(ids), error('No tissue voxels in selected ROIs.'); end
% Actual existing QC configuration, without saving any production files.
[~,~,~,~,~,~,oldBad,Q]=computeCESTQualityMetrics(d,...
    'RoiFile',roiFiles,'SaveTxt',false,'SaveMat',false,'SaveMaps',false);
th=Q.thresholds;
if th.UseRMSEPercentile, error('Absolute RMSE QC required.'); end
n=numel(ids);
Row=zeros(n,1); Col=Row; G=nan(n,1); Gcheck=G;
RM=nan(n,3); R2=RM; Zlo=RM; Zhi=RM;
Valid=false(n,3); Pass=Valid;
Bound=false(n,1); Spatial=Bound; ExistingBad=Bound;
Frac=nan(n,2); AbsEnd=Frac;
WorstPPM=nan(n,1);
for j=1:n
    [r,c]=ind2sub([nx ny],ids(j)); Row(j)=r; Col(j)=c;
    G(j)=F.gluc_integral_fullMap(r,c);
    ExistingBad(j)=oldBad(r,c);
    Bound(j)=Q.masks.Bound(r,c); Spatial(j)=Q.masks.Spatial(r,c);
    par=double(reshape(F.popt(r,c,1,:),[],1));
    if any(~isfinite(par)) || ~isfinite(G(j)), continue; end
    Gcheck(j)=par(17)*par(18)*pi/2;
    y=double(reshape(F.Z_corrExt(r,c,1,:),[],1));
    [fun,~]=fitmodelfunc_ANA(par,F.P);
    yh=double(reshape(fun(F.P.SEQ.w),[],1));
    if numel(yh)~=N, error('Model length mismatch at (%d,%d).',r,c); end
    residual=y-yh;
    if all(isfinite(residual))
        sse=sum(residual.^2);
        if sse>0, Frac(j,:)=[residual(kn)^2 residual(kp)^2]/sse; end
        AbsEnd(j,:)=[abs(residual(kn)) abs(residual(kp))];
        [~,iw]=max(abs(residual)); WorstPPM(j)=W(iw);
    end
    for v=1:3
        keep=true(N,1);
        if v==2, keep(kn)=false; elseif v==3, keep(kp)=false; end
        a=y(keep); b=yh(keep);
        Valid(j,v)=all(isfinite(a))&&all(isfinite(b));
        if ~Valid(j,v), continue; end
        RM(j,v)=sqrt(mean((a-b).^2));
        sst=sum((a-mean(a)).^2);
        if sst>0, R2(j,v)=1-sum((a-b).^2)/sst; end
        Zlo(j,v)=min(a); Zhi(j,v)=max(a);
        Pass(j,v)=isfinite(R2(j,v)) && ...
            RM(j,v)<=th.RMSEThreshUsed && R2(j,v)>=th.R2Thresh && ...
            Zlo(j,v)<=th.ZminThresh && Zhi(j,v)<=th.ZmaxThresh && ...
            ~Bound(j)&&~Spatial(j);
    end
end
% Rescue candidates require an originally evaluable full spectrum.
OriginalValidPass=~ExistingBad & Valid(:,1) & isfinite(G) & isfinite(Gcheck) & isfinite(R2(:,1));
CandidateNeg=ExistingBad & Valid(:,1) & Pass(:,2);
CandidatePos=ExistingBad & Valid(:,1) & Pass(:,3);
Category=repmat({'neither'},n,1);
ix=CandidateNeg & ~CandidatePos; Category(ix)=repmat({'negative_only'},nnz(ix),1);
ix=CandidatePos & ~CandidateNeg; Category(ix)=repmat({'positive_only'},nnz(ix),1);
ix=CandidateNeg & CandidatePos; Category(ix)=repmat({'either_endpoint'},nnz(ix),1);
Category(OriginalValidPass)=repmat({'already_pass'},nnz(OriginalValidPass),1);
ix=~Valid(:,1)|~isfinite(R2(:,1)); Category(ix)=repmat({'invalid_or_undefined_full_fit'},nnz(ix),1);
SessionPath=repmat({d},n,1);
T=table(SessionPath,Row,Col,G,Gcheck,ExistingBad,OriginalValidPass,Bound,Spatial,...
    repmat(negPPM,n,1),repmat(posPPM,n,1),WorstPPM,...
    'VariableNames',{'SessionPath','Row','Col','G_saved','G_from_parameters',...
    'ExistingBad','OriginalValidPass','BoundFlag','SpatialFlag',...
    'NegativeEndpointPPM','PositiveEndpointPPM','WorstResidualPPM'});
names={'Full','DropNeg','DropPos'};
for v=1:3
    nm=names{v};
    T.(['Valid_' nm])=Valid(:,v); T.(['RMSE_' nm])=RM(:,v);
    T.(['R2_' nm])=R2(:,v); T.(['Zmin_' nm])=Zlo(:,v);
    T.(['Zmax_' nm])=Zhi(:,v); T.(['Pass_' nm])=Pass(:,v);
end
T.NegativeSSEFraction=Frac(:,1); T.PositiveSSEFraction=Frac(:,2);
T.AbsNegativeResidual=AbsEnd(:,1); T.AbsPositiveResidual=AbsEnd(:,2);
T.CandidateNeg=CandidateNeg; T.CandidatePos=CandidatePos;
T.CandidateEither=CandidateNeg|CandidatePos; T.Category=Category;
U=table();
for k=1:numel(roiFiles)
    m=masks{k}; selected=m(ids);
    [~,label]=fileparts(roiFiles{k});
    u=table({d},{label},nnz(m),nnz(m&tissue),nnz(selected & Valid(:,1)),...
        nnz(selected & ExistingBad),nnz(selected & OriginalValidPass),...
        nnz(selected & CandidateNeg & ~CandidatePos),...
        nnz(selected & CandidatePos & ~CandidateNeg),...
        nnz(selected & CandidateNeg & CandidatePos),...
        nnz(selected & (CandidateNeg|CandidatePos)),...
        negPPM,posPPM,th.RMSEThreshUsed,th.R2Thresh,th.ZminThresh,th.ZmaxThresh,...
        'VariableNames',{'SessionPath','ROI','N_ROI','N_ROI_Tissue','N_ValidFull',...
        'N_ExistingBad','N_OriginalValidPass','N_NegativeOnly','N_PositiveOnly',...
        'N_EitherEndpoint','N_CandidateUnion','NegativeEndpointPPM','PositiveEndpointPPM',...
        'RMSEThreshold','R2Threshold','ZminThreshold','ZmaxThreshold'});
    if isempty(U), U=u; else, U=[U;u]; end %#ok<AGROW>
end
end

function sessions=local_sessions(root)
if ~exist(root,'dir'), error('Directory not found: %s',root); end
if exist(fullfile(root,'LorentzianFitResult.mat'),'file')
    sessions={root}; return;
end
D=[dir(fullfile(root,'*CEST_base*'));dir(fullfile(root,'*CEST_infu*'))];
D=D([D.isdir]);
sessions={};
for j=1:numel(D), sessions{end+1,1}=fullfile(D(j).folder,D(j).name); end %#ok<AGROW>
if isempty(sessions)
    S=[dir(fullfile(root,'*MIBIO*'));dir(fullfile(root,'*KHTCH*'))];
    S=S([S.isdir]);
    for k=1:numel(S)
        d=fullfile(S(k).folder,S(k).name);
        D=[dir(fullfile(d,'*CEST_base*'));dir(fullfile(d,'*CEST_infu*'))];
        D=D([D.isdir]);
        for j=1:numel(D), sessions{end+1,1}=fullfile(D(j).folder,D(j).name); end %#ok<AGROW>
    end
end
sessions=unique(sessions,'stable');
end
