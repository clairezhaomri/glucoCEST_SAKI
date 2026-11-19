# CEST QC — synced version (MATLAB R2018a)

## Files
| File | Status | What changed |
|---|---|---|
| `CESTQCConfig.m` | NEW | The only place where QC thresholds live. |
| `computeCESTQualityMetrics.m` | UPDATED | Defaults read from `CESTQCConfig`; thresholds + flag counts + per-ROI detail written to `QC_summary.txt`; also writes `QC_result.mat`; new optional `SpatialWindow`, `SaveMat`. Output signature and all variable names unchanged. |
| `plotCESTQCBeforeAfter.m` | UPDATED | No thresholds of its own (never had any); `RoiFile` default from config; threshold string printed at the bottom of the figure; warns if `BadPixelMask` given without `QCFlags`; interactive call no longer overwrites `QC_summary.txt`/`QC_result.mat`. |
| `WZ_QueueCEST_QC.m` | UPDATED | `roiFiles` from config; `forceRerun` / `makePlots` switches; passes `qcFlags` to the plot (same masks -> same reasons); saves `CESTQCConfig_used_<stamp>.txt` next to the CSV; CSV gains raw ROI mean/SD and pixel counts. |
| `WZ_QueueCEST_Analyze.m` | UPDATED (tiny) | `roiFiles` from config; fixed `fileparts` label bug on the error path. |
| `runCESTGlucoseROIFit.m (unchanged: QC is computed FROM its output, so it cannot apply QC itself)` | unchanged | — |
| `plotCESTGlucoseROIResults.m` | unchanged | — |

Install: copy all five changed/new files into `/Users/clairezhao/Documents/MATLAB/CEST/`, then in MATLAB:
```matlab
clear functions
which CESTQCConfig computeCESTQualityMetrics plotCESTQCBeforeAfter
```

## Where thresholds are set
Only in `CESTQCConfig.m`. Both the batch queue and interactive plotting read it.

Interactive experiment (does NOT touch saved QC files):
```matlab
plotCESTQCBeforeAfter(dataDir, 'QCOptions', {'ZminThresh', 0.2, 'R2Thresh', 0.7});
```
When you are happy, write the numbers into `CESTQCConfig.m`, set `forceRerun = true` in
`WZ_QueueCEST_QC.m`, run it, set `forceRerun = false` again.

## Exact names (do not rename)
Parameters (name-value AND fields of `CESTQCConfig`):
`RoiFile, ZminThresh, ZmaxThresh, UseRMSEPercentile, RMSEPercentile, RMSEThresh,
R2Thresh, UseBoundCheck, BoundTol, UseSpatialOutlier, SpatialMADMultiplier, SpatialWindow`

Outputs of `computeCESTQualityMetrics`:
`meanRMSE, meanR2, badFraction, roi_mean_clean, roi_std_clean, roi_labels, badPixelMask, qcFlags`

Internal maps: `Zmin_map, Zmax_map, RMSE_map, R2_map, atBound_map, tissueMask`

`qcFlags` fields: `masks.{Zmin,Zmax,RMSE,R2,Bound,Spatial}`, `values.{Zmin,Zmax,RMSE,R2}`,
`thresholds.<parameter names> + RMSEThreshUsed`, `thresholdString`, `overrides`, `tissue`, `roi(k)`

## QC_summary.txt layout
```
CEST fit quality summary
Folder: ...
Generated: ...
FitResultFile: .../LorentzianFitResult.mat (modified ...)
QCFunction: computeCESTQualityMetrics | ConfigFile: .../CESTQCConfig.m

Mean RMSE (tissue): ...          <- legacy lines, unchanged format
Mean R2 (tissue): ...
Bad pixel fraction (tissue): ...

ROI: ctx | Mean: ... | SD: ...
ROI: medulla | Mean: ... | SD: ...

[QC criteria]
ZminThresh: 0.1            (flag if Zmin_map > ZminThresh)
ZmaxThresh: 1.2            (flag if Zmax_map > ZmaxThresh)
UseRMSEPercentile: 1
RMSEPercentile: 70
RMSEThresh: Inf
RMSEThreshUsed: 0.0123     (flag if RMSE_map > RMSEThreshUsed)
R2Thresh: 0.8              (flag if R2_map < R2Thresh)
UseBoundCheck: 1
BoundTol: 1e-13
UseSpatialOutlier: 1
SpatialMADMultiplier: 5
SpatialWindow: [3 3]
RoiFile: ROI_ctx.mat, ROI_medulla.mat
Overrides vs CESTQCConfig: none

[Tissue flags]
nTissue: 9455
nFitted: ...
Zmin: n (x%)  Zmax: ...  RMSE: ...  R2: ...  Bound: ...  Spatial: ...
nBad (union): n (x%)

[ROI detail]
ROI: ctx | nROI: 201 | nFitted: .. | nBad: .. | nClean: .. | MeanRaw: .. | SDRaw: .. | MeanClean: .. | SDClean: ..
ROI: ctx flags | Zmin: .. | Zmax: .. | RMSE: .. | R2: .. | Bound: .. | Spatial: ..
```

## Python / Jupyter
`QC_result.mat` is saved with `-v7`:
```python
from scipy.io import loadmat
r = loadmat('QC_result.mat', squeeze_me=True, struct_as_record=False)
bad = r['badPixelMask']; thr = r['qcFlags'].thresholds; masks = r['qcFlags'].masks
```


## Post-QC mean-ROI Z-spectrum fit (plotCESTGlucoseROIResults.m, UPDATED)

Previously the mean ROI Z-spectrum was averaged over ALL ROI pixels, including
the ones QC rejects. Now:

| Situation | Behaviour |
|---|---|
| `plotCESTGlucoseROIResults(dataDir)` before any QC | pre-QC only, exactly as before |
| `plotCESTGlucoseROIResults(dataDir)` after QC ran | auto-loads `dataDir/QC_result.mat` (`'UseSavedQC'`, default true) |
| `... , 'BadPixelMask', badPixelMask, 'QCFlags', qcFlags` | uses exactly the masks you pass (this is what `WZ_QueueCEST_QC.m` does) |

When a QC mask is available, per ROI:

* `meanZspec_pre`  = mean over `ROI`                  -> fit -> `popt_roi_pre`
* `meanZspec_post` = mean over `ROI & ~badPixelMask`  -> fit -> `popt_roi_post`
* `gluc_integral_roiFit_pre/post` = `popt_roi(17) * popt_roi(18) * pi/2`
  (same definition as `gluc_integral_fullMap`)
* figures: `<roi>_MeanROI_Zspectrum_with_LorentzFit.png` (pre, as before),
  `<roi>_MeanROI_Zspectrum_with_LorentzFit_postQC.png`,
  `<roi>_MeanROI_Zspectrum_preVsPostQC.png` (QC thresholds printed at the bottom)
* `dataDir/LorentzianFitPostQC.mat` (`-v7`): `roiFit` struct array, `badPixelMask`,
  `qcThresholds`, `thresholdString`, `qcSource`, `generated`

New 4th output `roiFit`; outputs 1-3 (`roi_mean, roi_std, roi_labels`) keep
their old (pre-QC) meaning so `WZ_QueueCEST_Analyze.m` is unaffected.

`WZ_QueueCEST_QC.m` calls it after `plotCESTQCBeforeAfter` with the SAME
`badPixelMask/qcFlags` (`fitPostQC = true`, `'PlotMaps', false`) and adds CSV
columns `ROIfit_gluc_integral_pre`, `ROIfit_gluc_integral_post`.

### Batch CSV layout (CESTQCResults_<stamp>.csv)
One row per `AnimalID x session x ROI`:
`AnimalID, group (MIBIO|KHTCH), session (base|infu), sessionFolder, ROI, meanRMSE, meanR2, badFraction, ROI_n, ROI_nClean, ROI_mean_raw, ROI_std_raw, ROI_mean_clean, ROI_std_clean, ROIfit_gluc_integral_pre, ROIfit_gluc_integral_post, stamp`

```python
df = pd.read_csv('CESTQCResults_<stamp>.csv')
wide = df.pivot_table(index=['AnimalID','group','ROI'], columns='session',
                      values=['ROI_mean_clean','ROIfit_gluc_integral_post'])
wide.columns = [f'{m}_{s}' for m, s in wide.columns]
wide['delta_clean'] = wide['ROI_mean_clean_infu'] - wide['ROI_mean_clean_base']
```

Python:
```python
import scipy.io as sio
F = sio.loadmat('LorentzianFitPostQC.mat', squeeze_me=True, struct_as_record=False)
for r in F['roiFit']:
    print(r.label, r.nPix_pre, r.nPix_post, r.gluc_integral_roiFit_pre, r.gluc_integral_roiFit_post)
```
