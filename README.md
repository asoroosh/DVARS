# Introduction
Series of codes (mainly in MATLAB and Shell scripts), available here, support methods proposed in __*insight and
inference for DVARS* by Afyouni & Nichols__. You can access the preprint via xxxx

The toolbox can be used to:

* generate DVARS and five stanrdised variants of the measure.
* generate p-values for spikes to facilitate the descion whether a spike is corrupted and should be scrubbed. 
* explain variance of 4D fMRI images via three variance components: fast (D-var), slow (S-var) and edge (E-var).
* explain what share of the whole variance each component occupies.

# Configuration

## Dependencies
* MATLAB2016b (or higher) statistical toolbox [required]. 
* `'/Nifti_Util'` for neuroimaging analysis [optional].
* `FSL 5.0.9` (or newer) to produce DSE variance images [optional].     

## DVARS
Using `DVARSCalc.m` allows you to estimate the DVARS related statistics (e.g. D-var, p-values for spikes and standardised variants of the fast variance). 

As a quick example, for a pure white noise with *I=1000* observations (voxels) and *T=100* data-points, we have:  
```
Y = randn(1000,100);
[DVARS,Stat] = DVARSCalc(Y);
```

```
-Input is a Matrix.
-Extra-cranial areas removed: 1000x100
-No normalisation/scaling has been set!
-Data centred. Untouched Grand Mean: 0.0019315, Post-norm Grand Mean: 0.0019315
-Settings: 
--Test Method:          Z
--ExpVal method:        median
--VarEst method:        hIQRd
--Power Transformation: 1

Settings: TestMethod=Z  I=1000  T=100 
----Expected Values----------------------------------
    sig2bar    sig2median    median    sigbar2     xbar 
    _______    __________    ______    _______    ______

    2.0342     1.984         1.9898    2.004      1.9998

----Variances----------------------------------------
     S2       IQRd         hIQRd  
    ____    _________    _________

    0.01    0.0052876    0.0034822
```

As you can see in the printed log of the function, the default method of testing is Z of expected value `median` and variance of `hIQR`. 
However, these settings are easily amendable. See example below:

```
[DVARS,Stat]=DVARSCalc(Y,'TestMeth','X2','VarRobIQR','hIQR','TransPower',1/3);
```

## DSE Variance Decomposition


## Visualisation


