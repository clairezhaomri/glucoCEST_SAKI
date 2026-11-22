function WZ_QueueCEST_ControlQC_PDF()
% Read existing PNGs only; no fitting, QC computation, or source modification.
% MATLAB R2021b or newer.
Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
ReportBase = fullfile(Basepath,'QC_Report');
SavePrefix = '';
ROILabels = {'ctx','medulla'};
GroupOrder = {'Ctrl+SC','Ctrl+HFD','Ctrl+HPD'};
% Only needed if there are repeated baseline/infusion sessions:
% {exact_subject_folder, exact_baseline_folder, exact_infusion_folder}
% Empty session folder name means missing side.
SessionOverride = cell(0,3);

assert(isfolder(Basepath),'Basepath not found: %s',Basepath);
assert(~verLessThan('matlab','9.11'),'Requires MATLAB R2021b or newer.');
if ~isfolder(ReportBase),mkdir(ReportBase);end
stamp=datestr(now,'yyyymmdd_HHMMSS');
RunDir=fullfile(ReportBase,['ControlQC_' stamp]);
k=1;
while isfolder(RunDir)
    RunDir=fullfile(ReportBase,sprintf('ControlQC_%s_%02d',stamp,k));k=k+1;
end
mkdir(RunDir);
PDF=fullfile(RunDir,'Controls_baseline_vs_infusion.pdf');
fprintf('Output directory:\n%s\n',RunDir);
D=[dir(fullfile(Basepath,'*MIBIO*'));dir(fullfile(Basepath,'*KHTCH*'))];
D=D([D.isdir]);
if ~isempty(D),[~,ix]=unique({D.name},'stable');D=D(ix);end
inventory=cell(0,7);selected=cell(0,3);
for i=1:numel(D)
    name=D(i).name;
    isCtrl=contains(name,'Ctrl','IgnoreCase',false);
    hits=[contains(name,'Normal','IgnoreCase',false), ...
          contains(name,'Fat','IgnoreCase',false), ...
          contains(name,'Protein','IgnoreCase',false)];
    group='';
    if ~isCtrl
        status='NOT_CONTROL';
    elseif sum(hits)==0
        status='CONTROL_DIET_UNRECOGNIZED';
    elseif sum(hits)>1
        status='AMBIGUOUS_DIET';
    else
        rank=find(hits,1);group=GroupOrder{rank};status='SELECTED';
        selected(end+1,:)={group,name,rank};
    end
    root=fullfile(Basepath,name);
    b=local_sessions(root,'*CEST_base*');u=local_sessions(root,'*CEST_infu*');
    inventory(end+1,:)={name,isCtrl,group,status,numel(b),numel(u),root};
end
I=cell2table(inventory,'VariableNames',{'SubjectFolder','ContainsCtrl','Group','SelectionStatus','NBaseline','NInfusion','SubjectPath'});
writetable(I,fullfile(RunDir,'subject_inventory.csv'));
if isempty(selected)
    warning('No controls selected. See subject_inventory.csv. No PDF created.');return
end
T=cell2table(selected,'VariableNames',{'Group','SubjectFolder','GroupRank'});
T=sortrows(T,{'GroupRank','SubjectFolder'});
writetable(T,fullfile(RunDir,'selected_controls.csv'));
manifest=cell(0,8);page=0;
for i=1:height(T)
    group=T.Group{i};subject=T.SubjectFolder{i};root=fullfile(Basepath,subject);
    b=local_sessions(root,'*CEST_base*');u=local_sessions(root,'*CEST_infu*');
    ov=find(strcmp(SessionOverride(:,1),subject));
    assert(numel(ov)<=1,'Duplicate session override for %s',subject);
    if ~isempty(ov)
        b=local_override(root,SessionOverride{ov,2});
        u=local_override(root,SessionOverride{ov,3});
    end
    if numel(b)>1 || numel(u)>1
        manifest(end+1,:)={group,subject,'',strjoin(b,' | '),strjoin(u,' | '),NaN,'AMBIGUOUS_SESSIONS','Set SessionOverride; no automatic pairing'};
        local_manifest(manifest,RunDir);
        warning('Skipped repeated sessions: %s. Specify SessionOverride.',subject);continue
    end
    bn='';un='';if ~isempty(b),bn=b{1};end;if ~isempty(u),un=u{1};end
    for r=1:numel(ROILabels)
        roi=ROILabels{r};name=[SavePrefix roi '_QC_before_after.png'];
        bp='';up='';
        if ~isempty(bn),bp=fullfile(root,bn,name);end
        if ~isempty(un),up=fullfile(root,un,name);end
        hasB=~isempty(bp)&&isfile(bp);hasU=~isempty(up)&&isfile(up);
        if ~hasB && ~hasU
            manifest(end+1,:)={group,subject,roi,bp,up,NaN,'BOTH_PNG_MISSING',''};
            local_manifest(manifest,RunDir);continue
        end
        try
            f=figure('Visible','off','Color','w','Units','pixels','Position',[50 50 2800 1100]);
            figCleanup=onCleanup(@()local_close(f));
            annotation(f,'textbox',[.01 .94 .98 .055], ...
                'String',sprintf('%s | %s | ROI: %s',group,subject,roi), ...
                'Interpreter','none','EdgeColor','none','HorizontalAlignment','center','FontSize',15);
            local_panel(f,[.012 .08 .478 .80],bp,hasB,['BASELINE | ' bn]);
            local_panel(f,[.510 .08 .478 .80],up,hasU,['AFTER GLUCOSE INFUSION | ' un]);
            annotation(f,'textbox',[.015 .005 .97 .05], ...
                'String','Original PNGs only. Thresholds, labels and clipping unchanged. Missing panel does not indicate QC failure.', ...
                'Interpreter','none','EdgeColor','none','FontSize',11);
            drawnow;
        catch ME
            if exist('figCleanup','var'),clear figCleanup;end
            manifest(end+1,:)={group,subject,roi,bp,up,NaN,'PNG_RENDER_ERROR',ME.message};
            local_manifest(manifest,RunDir);
            warning('Render error %s / %s: %s',subject,roi,ME.message);continue
        end
        try
            exportgraphics(f,PDF,'ContentType','image','Resolution',200,'Append',page>0);
        catch ME
            clear figCleanup
            manifest(end+1,:)={group,subject,roi,bp,up,NaN,'PDF_EXPORT_FAILED',ME.message};
            local_manifest(manifest,RunDir);
            error('PDF export failed; partial PDF may exist. See manifest. %s',ME.message);
        end
        page=page+1;clear figCleanup
        status='PAIRED';
        if ~hasB,status='BASELINE_PNG_MISSING';elseif ~hasU,status='INFUSION_PNG_MISSING';end
        manifest(end+1,:)={group,subject,roi,bp,up,page,status,''};
        local_manifest(manifest,RunDir);
        fprintf('Page %d: %s | %s | %s\n',page,group,subject,roi);
    end
end
local_manifest(manifest,RunDir);
for j=1:numel(GroupOrder)
    fprintf('%s: %d selected subject folders\n',GroupOrder{j},sum(strcmp(T.Group,GroupOrder{j})));
end
if page>0
    fprintf('\nDone: %d pages\n%s\n',page,PDF);
else
    warning('No PDF pages created. Inspect manifest.csv.');
end
end

function names=local_sessions(root,pattern)
d=dir(fullfile(root,pattern));d=d([d.isdir]);names=sort({d.name});
end
function names=local_override(root,name)
if isempty(name),names={};return;end
assert(isfolder(fullfile(root,name)),'Override session does not exist: %s',fullfile(root,name));
names={name};
end
function local_manifest(rows,root)
M=cell2table(rows,'VariableNames',{'Group','SubjectFolder','ROI','BaselinePNG','InfusionPNG','Page','Status','Message'});
writetable(M,fullfile(root,'manifest.csv'));
end
function local_close(f)
if isgraphics(f),close(f);end
end
function local_panel(f,pos,path,hasImage,label)
ax=axes('Parent',f,'Position',pos);
if hasImage
    [im,map]=imread(path);
    if ~isempty(map),im=ind2rgb(im,map);end
    image(ax,im);axis(ax,'image');axis(ax,'off');
else
    axis(ax,[0 1 0 1]);axis(ax,'off');
    text(ax,.5,.5,'MISSING PNG / SESSION','HorizontalAlignment','center','FontSize',18,'Color',[.7 0 0]);
end
title(ax,label,'Interpreter','none','FontSize',12,'FontWeight','normal');
end
