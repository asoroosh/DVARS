# DSE & DVARS (Python)

Python implementation of DSE variance decomposition and DVARS inference for fMRI data.

Reference: Afyouni, S. & Nichols, T.E. "Insight and Inference for DVARS." *NeuroImage* 172 (2018): 291-312.

## Requirements

```
numpy, scipy, nibabel, matplotlib
```

## Usage

```python
from DSE import DSE_Calc, DVARS_Calc, CleanNIFTI

# Y can be a path to a NIfTI file, or a 2D (voxels x time) / 4D array
DSEOut   = DSE_Calc(Y)
DVARSOut = DVARS_Calc(Y)
```

See `DSEnDVAR_Demo.py` for a worked example.

## Functions

- **`CleanNIFTI(V, scl=0, norm=0, demean=True, copy=True)`** — Loads/reshapes input (NIfTI path, or 2D/4D array) into an `I x T` matrix, drops zero/NaN voxels, and optionally scales/normalises and demeans. Returns `(Y, T, I, I1, rmvIdx, imobj, GS, MeanImage)`.
- **`DSE_Calc(Y, NewDir2DSE='', DSEImage=False, verbose=False, ...)`** — Decomposes variance into Averaged (A), Delta (D), Slow (S), and Edge (E) components, at whole, global, and non-global levels. Returns a dict with scalar variances, time series, percentages, and an ANOVA-style table. Set `DSEImage=True` to also write 3D/4D variance NIfTI images (requires the `imobj` from `CleanNIFTI`), or `verbose=True` to save tables/plots to `NewDir2DSE`.
- **`DVARS_Calc(Y, dd=1, WhichExpVal='median', WhichVar='hIQRd', DeltapDvarThr=5, ...)`** — Computes standardised DVARS (`NDVARS_X2`) and flags statistically/practically significant volumes via chi-squared inference. Returns a dict of `DVARS`, `Inference`, and `Stats`.
- **`SaveMe2txt` / `SaveMe2Nifti`** — Helpers to write result dicts to `.txt` or `.nii.gz` files.
