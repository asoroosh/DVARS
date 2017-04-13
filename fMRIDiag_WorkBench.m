clear

%replace this path with the DVARS directory on your machine
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/DVARS'))

%replace this path with the path to the nifti/cifti images on your machine
V0=['/Volumes/HCP_S900/HCP_10Unrel_Vols/115320/115320/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz'];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the simplest form of using DSE variance decomposition and DVARS
% inference on the data. As the example is for HCP data, we used 1/100
% downscaling factor. The scale/norm factors should be changed for other
% data-sets (However, if you have processed your data via FSL, leave the 
% scaling 1/100 as it is).

%--NIFTI--------------------------------------------
V1 = load_untouch_nii(V0);
V2 = V1.img; 
X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
I0 = prod([X0,Y0,Z0]);
Y  = reshape(V2,[I0,T0]); clear V2 V1; 
% %--CIFTI--------------------------------------------
% addpath /Users/sorooshafyouni/Home/matlab/Ext/cifti-matlab-master
% V1 = ft_read_cifti(V0);
% V2=V1.dtseries;
% I0=size(V2,1); T0=size(V2,2);
% Y=V2; clear V2 V1; 
% %

[DVARS,DVARS_Stat]=DVARSCalc(Y,'scale',1/100);
[V,DSE_Stat]=DSEvars(Y,'scale',1/100);

fMRIDiag_plot(V,DVARS_Stat)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this example, you can see how to add BOLD intensity images to your
% figure to have a better understanding the of the effect of spikes/excursions 
% in your variance components. Read help section of 'fMRIDiag_plot' for
% further parameters that can be used. 

%--NIFTI--------------------------------------------
V1 = load_untouch_nii(V0);
V2 = V1.img; 
X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
I0 = prod([X0,Y0,Z0]);
Y  = reshape(V2,[I0,T0]); clear V2 V1; 
% %--CIFTI--------------------------------------------
% addpath /Users/sorooshafyouni/Home/matlab/Ext/cifti-matlab-master
% V1 = ft_read_cifti(V0);
% V2=V1.dtseries;
% I0=size(V2,1); T0=size(V2,2);
% Y=V2; clear V2 V1; 
% %

[DVARS,DVARS_Stat]=DVARSCalc(Y,'scale',1/100);
[V,DSE_Stat]=DSEvars(Y,'scale',1/100);

fMRIDiag_plot(V,DVARS_Stat,'BOLD',Y,'ColRng',[0 100])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this example, it is assumed that the desplacement information has also
% been provided for a subject. An extra sub-plot, automatically, is added
% to the top of your figure to report the FD and Absolute values of the
% movements acquired via 'FDCalc'. 

%--NIFTI--------------------------------------------
V1 = load_untouch_nii(V0);
V2 = V1.img; 
X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
I0 = prod([X0,Y0,Z0]);
Y  = reshape(V2,[I0,T0]); clear V2 V1; 
% %--CIFTI--------------------------------------------
% addpath /Users/sorooshafyouni/Home/matlab/Ext/cifti-matlab-master
% V1 = ft_read_cifti(V0);
% V2=V1.dtseries;
% I0=size(V2,1); T0=size(V2,2);
% Y=V2; clear V2 V1; 
% %

%--Movement Parameters--------------------------------
%replace this path with the path to the text files with movement regressors
%of the image on your machine
%Note that MovPartextImport only works safe with the HCP files, you have to
%insert the movement regressors mannually (just drag the text file into the
%workspace!) to the Matlab.

MovPar=MovPartextImport(['/Volumes/HCP_S900/HCP_10Unrel_Vols/115320/115320/MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt']);
[FDts,FD_Stat]=FDCalc(MovPar);
%

[DVARS,DVARS_Stat]=DVARSCalc(Y,'scale',1/100);
[V,DSE_Stat]=DSEvars(Y,'scale',1/100);

fMRIDiag_plot(V,DVARS_Stat,'BOLD',Y,'FD',FDts,'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans])
