%% Archive old (pre-correction) fit results before the full rerun.
% Moves LorentzianFitResult.mat, QC_summary.txt, and *_QC_before_after.png
% out of every subject folder into a single archive directory, so
% WZ_QueueCEST_Analyze.m's existing "skip if already fit" check naturally
% triggers a refit for every subject using the corrected data -- without
% editing the skip logic itself, and without losing the old results in
% case you want to compare before/after correction later.

Basepath = '/Users/clairezhao/Documents/data/Hoerr/MRI/pv6/CEST/nutritionCEST/NewData/';
ArchiveDir = fullfile(Basepath, '_archived_precorrection_fits');
if ~isfolder(ArchiveDir), mkdir(ArchiveDir); end

filesToArchive = dir(fullfile(Basepath, '**', 'LorentzianFitResult.mat'));
filesToArchive = [filesToArchive; dir(fullfile(Basepath, '**', 'QC_summary.txt'))];
filesToArchive = [filesToArchive; dir(fullfile(Basepath, '**', '*_QC_before_after.png'))];

for k = 1:numel(filesToArchive)
    srcPath = fullfile(filesToArchive(k).folder, filesToArchive(k).name);
    [~, subjFolderName] = fileparts(filesToArchive(k).folder);
    destSub = fullfile(ArchiveDir, subjFolderName);
    if ~isfolder(destSub), mkdir(destSub); end
    movefile(srcPath, fullfile(destSub, filesToArchive(k).name));
end

fprintf('Archived %d files to %s\n', numel(filesToArchive), ArchiveDir);
