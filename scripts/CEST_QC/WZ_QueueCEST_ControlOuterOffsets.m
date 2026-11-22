function [All, Summary, runDir] = WZ_QueueCEST_ControlOuterOffsets(Basepath)
% Controls-only diagnostic. Requires testCESTOuterOffsetsBatch.m.
% Baseline: process if exactly one *CEST_base* directory exists.
% Infusion: process ONLY the exact case-sensitive directory CEST_infu.
% Call the core separately on explicit session paths to prevent wildcard leakage.
if nargin<1 || isempty(Basepath)
    Basepath='/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
end
assert(isfolder(Basepath),'Basepath not found: %s',Basepath);
assert(exist('testCESTOuterOffsetsBatch','file')~=0,'Missing testCESTOuterOffsetsBatch.m');
reportBase=fullfile(Basepath,'QC_Report');
if ~isfolder(reportBase), mkdir(reportBase); end
[~,token]=fileparts(tempname(reportBase));
runDir=fullfile(reportBase,['ControlOuterOffsets_' datestr(now,'yyyymmdd_HHMMSS') '_' token]);
mkdir(runDir);
groups={'Ctrl+SC','Ctrl+HFD','Ctrl+HPD'};
D=[dir(fullfile(Basepath,'*MIBIO*'));dir(fullfile(Basepath,'*KHTCH*'))];
D=D([D.isdir]);
if ~isempty(D), [~,ix]=unique({D.name},'stable'); D=D(ix); end
inventory=cell(0,7); selected=cell(0,3);
for i=1:numel(D)
    name=D(i).name; root=fullfile(Basepath,name);
    ctrl=contains(name,'Ctrl','IgnoreCase',false);
    hits=[contains(name,'Normal','IgnoreCase',false),...
        contains(name,'Fat','IgnoreCase',false),contains(name,'Protein','IgnoreCase',false)];
    group='';
    if ~ctrl, status='NOT_CONTROL';
    elseif sum(hits)==0, status='CONTROL_DIET_UNRECOGNIZED';
    elseif sum(hits)>1, status='AMBIGUOUS_DIET';
    else
        rank=find(hits,1); group=groups{rank}; status='SELECTED';
        selected(end+1,:)={group,name,rank}; %#ok<AGROW>
    end
    b=local_dirs(root,'*CEST_base*'); u=local_dirs(root,'*CEST_infu*');
    inventory(end+1,:)={name,ctrl,group,status,numel(b),numel(u),root}; %#ok<AGROW>
end
writetable(cell2table(inventory,'VariableNames',{'SubjectFolder','ContainsCtrl',...
    'Group','SelectionStatus','NBaseline','NInfusion','SubjectPath'}),...
    fullfile(runDir,'subject_inventory.csv'));
S=cell2table(selected,'VariableNames',{'Group','SubjectFolder','GroupRank'});
if ~isempty(S), S=sortrows(S,{'GroupRank','SubjectFolder'}); end
writetable(S,fullfile(runDir,'selected_controls.csv'));
All=table(); Summary=table(); manifest=cell(0,7);
for i=1:height(S)
    name=S.SubjectFolder{i}; group=S.Group{i}; root=fullfile(Basepath,name);
    b=local_dirs(root,'*CEST_base*'); u=local_dirs(root,'*CEST_infu*');
    paths={}; phases={};
    if numel(b)==1
        paths{end+1}=fullfile(root,b{1}); phases{end+1}='base'; %#ok<AGROW>
    elseif isempty(b)
        manifest(end+1,:)={group,name,'base','','MISSING_BASELINE','',''}; %#ok<AGROW>
    else
        manifest(end+1,:)={group,name,'base',root,'AMBIGUOUS_BASELINE','',strjoin(b,' | ')}; %#ok<AGROW>
    end
    exact=find(strcmp(u,'CEST_infu'));
    if numel(exact)==1
        paths{end+1}=fullfile(root,u{exact}); phases{end+1}='infu'; %#ok<AGROW>
    else
        manifest(end+1,:)={group,name,'infu',fullfile(root,'CEST_infu'),...
            'MISSING_EXACT_CEST_INFU','','No substitute infusion folder selected.'}; %#ok<AGROW>
    end
    ignored=u(~strcmp(u,'CEST_infu'));
    for k=1:numel(ignored)
        manifest(end+1,:)={group,name,'infu',fullfile(root,ignored{k}),...
            'IGNORED_NONEXACT_INFUSION','','Only exact CEST_infu is eligible.'}; %#ok<AGROW>
    end
    for k=1:numel(paths)
        d=paths{k}; phase=phases{k};
        fprintf('\nControl %d/%d: %s | %s | %s\n',i,height(S),group,name,phase);
        if ~isfile(fullfile(d,'LorentzianFitResult.mat'))
            manifest(end+1,:)={group,name,phase,d,'MISSING_FIT','','No refit performed.'}; %#ok<AGROW>
            continue;
        end
        try
            [A,B,child]=testCESTOuterOffsetsBatch(d,'OutputRoot',runDir);
            if ~isempty(A)
                A.Group=repmat({group},height(A),1);
                A.SubjectFolder=repmat({name},height(A),1);
                A.SessionType=repmat({phase},height(A),1);
                if isempty(All), All=A; else, All=[All;A]; end %#ok<AGROW>
            end
            if ~isempty(B)
                B.Group=repmat({group},height(B),1);
                B.SubjectFolder=repmat({name},height(B),1);
                B.SessionType=repmat({phase},height(B),1);
                if isempty(Summary), Summary=B; else, Summary=[Summary;B]; end %#ok<AGROW>
            end
            if isempty(A), status='SESSION_FAILED'; msg='See child session_log.csv.';
            else, status='OK'; msg=''; end
            manifest(end+1,:)={group,name,phase,d,status,child,msg}; %#ok<AGROW>
        catch ME
            manifest(end+1,:)={group,name,phase,d,'ERROR','',ME.message}; %#ok<AGROW>
            warning('ControlOuterOffsets:sessionFailed','%s: %s',d,ME.message);
        end
        local_save(All,Summary,manifest,runDir);
    end
    local_save(All,Summary,manifest,runDir);
end
local_save(All,Summary,manifest,runDir);
fprintf('\nDone. Only exact CEST_infu included; source fits/QC unchanged.\n%s\n',runDir);
end

function names=local_dirs(root,pattern)
d=dir(fullfile(root,pattern)); d=d([d.isdir]); names=sort({d.name});
end
function local_save(All,Summary,rows,root)
writetable(cell2table(rows,'VariableNames',{'Group','SubjectFolder','SessionType',...
    'SessionPath','Status','OutputDirectory','Message'}),fullfile(root,'control_manifest.csv'));
if ~isempty(All), writetable(All,fullfile(root,'CONTROLS_ALL_voxels.csv')); end
if ~isempty(Summary), writetable(Summary,fullfile(root,'CONTROLS_ALL_roi_summary.csv')); end
end
