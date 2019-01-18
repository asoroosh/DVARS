clear

V = '/Users/sorooshafyouni/Desktop/HCP/135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz';

Y = CleanNIFTI(V,'scale',1/100,'demean'); 

[DSE,DSE_Stats] = DSEvars(Y);

[DVARS,DVARS_Stats] = DVARSCalc(Y);

