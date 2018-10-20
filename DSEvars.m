function [V,Stat]=DSEvars(V0,varargin)
%[V,Stat]=DSEvars(V0,varargin)
%
%%%INPUTS:
%   V0:             Can be (1) a string indicating the path to the 
%                   nifti/cifti file (2) a numerical matrix of size IxT. 
%                   Where I is number of voxels (I=Nx x Ny x Nz) and T is 
%                   number of data-points.
% OPTIONS:
%
%
%   'DestDir'     : Output directory. Should only be used when the
%                   input is a nifti image and user needs to save the
%                   S, D and E (3D and 4D) images. 
%                   e.g.: [V,Stat]=DSEvars(V0,'DestDir','~/Where/to/save/')
%
%   'saveDSEtable': If triggered and followed by a path + filename.csv it
%                   saves the DSE table as a csv file
%
%   'Norm'        : Intensity normalisation to a given scale.
%                   e.g.: [V,Stat]=DSEvars(V0,'Norm',100)
%
%   'Scale'       : Scale the intensity between the data-sets.
%                   e.g.: [V,Stat]=DSEvars(V0,'Scale',1/10)
%
%   'verbose'     : Set to 1 if you need the log of runing code 
%                   [default:1]
%                   e.g.: [V,Stat]=DSEvars(V0,'verbose',1)
%%%OUTPUTS:
%
%   V:       Structure contains the variance components:
%               V.{A,S,D,E}var:    time series of var components
%               V.w_{A,S,D,E}var:  sum of mean squared of components
%               V.g_{A,S,D,E}var:  sum of mean squared of global components
%               V.ng_{A,S,D,E}var: sum of mean squared of non-global components
%   Stat:    Structure contains the higher level parameters of the comps:
%               Stat.Labels:    Labels indicating the order of next vars
%               Stat.SS:        Sum-squared 
%               Stat.MS:        Mean-squared
%               Stat.RMS:       Root-Mean-Squared
%               Stat.Prntg:     Percentage of the whole variance
%               Stat.RelVar:    Percentage of the whole variance relative
%                               to the iid case.
%
%               Stat.DeltapDvar: \Delta\%D-var
%               Stat.pDvar:      \%D-var
%               Stat.DeltapSvar: \Delta\%S-var
%               Stat.pSvar:      \%S-var
%
%
%%%NOTES:
%   1) It is recommended to only use time series of intra-cranial voxels as
%   inclusding the extra-cranial may inflate the variance. You can use
%   'bet' in FSL package to remove the extra-cranial areas. The scripts
%   automatically remove the zero/NaN voxels.
%
%   2) If a destination directory doesn't exist, the function automatically
%   make a directory with the given 'DestDir'.
%
%   3) To fully exploit the DSEvars, the data should *NOT* be undergone any
%   form of temporal filtering, as temporal filtering may remove the high
%   freq fluctuations.
%
%   4) For inter-site/cohort comparison, it is recommended that the
%   intensity is scale accordingly by option 'Norm' or 'Scale'.
%
%   5) If the input is set to be a NIFTI file, you require Nifti_Util 
%      (provided in the directory). For input of CIFTI you require to
%      addpath the FieldTrip toolbox from: 
%      http://www.fieldtriptoolbox.org/reference/ft_read_cifti 
%
%%%EXAMPLE:
%
%   For iid case, numerical matrix:
%
%   I=4e4; T=1200; Y=randn(I,T);
%   [V,Stat]=DSEvars(Y);
%   In this example, the function returns the variance components and print
%   the SS and ANOVA tables for input of numerical matrix.
%   
%   OneSub='~/100307/rfMRI_REST1_LR.nii.gz' %a HCP Subject
%   [V,Stat]=DSEvars(OneSub);
%   In this example, the function returns the variance components and print
%   the SS and ANOVA tables for input of nifti image. 
%   
%   OneSub='~/100307/rfMRI_REST1_LR.nii.gz' %a HCP Subject
%   [V,Stat]=DSEvars(OneSub,'verbose',1,'DestDir','~/temp','Norm',100);
%
%   Stat.DpDVARS  : \Delta\%D-var (Exceed fast Standardised DVARS) 
%   Stat.pDvar    : \%D-var       (Percentage of the whole var -A-var-) 
%   In this example, the function returns the variance components and print
%   the SS and ANOVA tables for input of nifti image. It also saves the 4D 
%   and 3D images of variance components in directory '~/temp'.     
%
%
%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021.1
%
%
%%%
%   Soroosh Afyouni & Thomas Nichols, UoW, Feb 2017
% 
%   https://github.com/asoroosh/DVARS
%   http://warwick.ac.uk/tenichols
%
%   Please report bugs to srafyouni@gmail.com
%_________________________________________________________________________
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________


%% ParCheck
t3_varn = {'Avar','Dvar','Svar','Evar'};
t3_rown = {'Whole','Global','non-Global'};
Row_labs = {'Avar','Dvar','Svar','Evar','g_Avar','g_Dvar','g_Svar','g_Evar'};
Col_labs = {'MS','RMS','Percentage_of_whole','Relative_to_iid'};

% Input Check-------------------------

gsrflag=0; verbose=1; DestDir=[]; DestDirTable=[]; md=[]; scl=[];
if sum(strcmpi(varargin,'gsrflag'))
   gsrflag      =   varargin{find(strcmpi(varargin,'gsrflag'))+1};
end
if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
end

if sum(strcmpi(varargin,'saveDSEtable'))
   DestDirTable      =   varargin{find(strcmpi(varargin,'saveDSEtable'))+1};
end

if sum(strcmpi(varargin,'destdir'))
   DestDir      =   varargin{find(strcmpi(varargin,'destdir'))+1};
   if sum(strcmpi(varargin,'images'))
      imagelist  =   varargin{find(strcmpi(varargin,'images'))+1}; 
   else
      imagelist  =   {'Dvar','Svar'};
   end
   %add something here to show the var images, just in case; with a verbose
   %trigger, of course!
end
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end

% Add toolbox to open the images-------
if isempty(strfind(path,'Nifti_Util'))
    if verbose; disp('-Nifti_Util added to the path.'); end;
    addpath(genpath('Nifti_Util'));
end

%---temp
% if sum(strcmpi(varargin,'MeanImage'))
%    mYr    =   varargin{find(strcmpi(varargin,'MeanImage'))+1};
%    mYr=mYr(mYr~=0 & ~isnan(mYr));
%    %size(mYr)
%    md           =   median(mYr);
% end

if ischar(V0)
    [ffpathstr,ffname,ffext]=fileparts(V0);
    if verbose; disp(['-Path to the image is: ' ffpathstr]); end;
    
    
    %if you are using MATLAB <2016, please replace 'contains' with 'strfind'
    if contains(ffname,'.dtseries') || contains(ffext,'.dtseries')
        if verbose; disp(['--File is CIFTI: ' ffname ffext]); end;
        V1=ft_read_cifti(V0);
        V2=V1.dtseries;
        I0=size(V2,1); T0=size(V2,2);
        Y=V2; clear V2 V1; 
    %if you are using MATLAB <2016, please replace 'contains' with 'strfind'    
    elseif ~contains(ffname,'.dtseries') || contains(ffname,'.nii') 
        if verbose; disp(['--File is NIFTI: ' ffname ffext]); end;
        V1 = load_untouch_nii(V0);
        V2 = V1.img; 
        X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
        I0 = prod([X0,Y0,Z0]);
        Y  = reshape(V2,[I0,T0]); clear V2;
    else
        error('Unknown input image.')
    end
    
    if verbose; disp('-Image loaded.'); end;
elseif isnumeric(V0) %&& size(V0,1)>size(V0,2)
    if verbose; disp('-Input is a Matrix.'); end;
    if size(V0,1)<=size(V0,2)
        warning('Check the input, matrix should be in form of IxT, where I=XxYxZ!');  
    end
    Y = double(V0);  
    I0= size(Y,1); T0 = size(Y,2);
%elseif isnumeric(V0) && size(V0,1)<=size(V0,2)
%    if verbose; disp('-Input is a Matrix.'); end;
%    warning('Check the input, matrix should be in form of IxT, where I=XxYxZ!');    
end

Y = double(Y);%to work with int 16bit as well.

mvY_WholeImage = mean(Y,2);
%Remove voxels of zeros/NaNs---------------------------------------------------
nan_idx    = find(isnan(sum(Y,2)));
zeros_idx  = find(sum(Y,2)==0);
idx        = 1:I0;
idx([nan_idx;zeros_idx]) = [];
Y([nan_idx;zeros_idx],:) = [];
I1 = size(Y,1); %update number of voxels
if verbose; disp(['-Extra-cranial areas removed: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]); end;

mvY_Untouched = mean(Y,2);

% Intensity Normalisation------------------------------------------------------
IntnstyScl = @(Y,md,scl) (Y./md).*scl; 
if ~isempty(scl) && isempty(md)
    md  = median(mean(Y,2)); %NB median of the mean image.
    %md  = mean(mean(Y,2)); %NB *mean* of the mean image.
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Normalised by ' num2str(scl) '&' num2str(md) '.']); end;
elseif ~isempty(scl) && ~isempty(md)
    assert(md==1,'4D mean in scalling cannot be anything other than 1!')
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Scaled by ' num2str(scl) '.']); end;
elseif isempty(scl) && isempty(md)    
    if verbose; disp('-No normalisation/scaling has been set!'); end;
else
    error('Something is wrong with param re: intensity normalisation')
end
%Centre the data-----------------------------
mvY_NormInt    =    mean(Y,2); %later will be used as grand mean! don't touch it!
dmeaner =    repmat(mvY_NormInt,[1,T0]);
Y       =    Y-dmeaner; clear dmeaner

mvY_Demeaned = mean(Y,2);
%----------
if verbose; disp(['-Data centred. Untouched Grand Mean: ' num2str(mean(mvY_Untouched)) ', Post-norm Grand Mean: ' num2str(mean(mvY_NormInt)) ', Post demean: ' num2str(mean(mvY_Demeaned))]); end;
%Data GSRed--------------------------------ONLY FOR TEST-----------------
if gsrflag 
    Y  =    fcn_GSR(Y);
    if verbose; disp('-Data GSRed.'); end;
end
%------------------------------------------ONLY FOR TEST-----------------
%% Lagged Images
B.Ybar  = sum(Y)./I1; %global signal is here!

D       = Y(:,1:end-1)-Y(:,2:end);
B.Dbar  = sum(D)./I1;

S       = Y(:,1:end-1)+Y(:,2:end);
B.Sbar  = sum(S)./I1;

Ytail   = Y(:,end); Yhead=Y(:,1);
B.Ytbar = sum(Ytail)./I1; 
B.Y1bar = sum(Yhead)./I1;
%% DSE Var Images
%4D images
V_Img.Avar_ts = Y.^2;
V_Img.Dvar_ts = D.^2./4;
V_Img.Svar_ts = S.^2./4;
V_Img.Evar_ts = [Yhead,Ytail].^2./2;
%3D images -- averaged across time. 
V_Img.Avar = mean(Y.^2,2);
V_Img.Dvar = mean(D.^2,2)./2;
V_Img.Svar = mean(S.^2,2)./2;
V_Img.Evar = mean([Yhead.^2,Ytail.^2],2); % <<<< should be checked
%% DSE Time series -- averaged across I
V.Avar_ts = mean(V_Img.Avar_ts);
V.Dvar_ts = mean(V_Img.Dvar_ts);
V.Svar_ts = mean(V_Img.Svar_ts);
V.Evar_ts = mean(V_Img.Evar_ts);
%% Save Images?
if ~isempty(DestDir) && ischar(V0)
    %if ~any(strfind(path,'spm')); warning('**SPM has not been added to the path!**'); end;
    if exist(DestDir,'dir')~=7; mkdir(DestDir); end;
    
    %savedir = [pwd '/' DestDir '/']; 
    for is=imagelist
        if verbose; disp(['****' is{1} ':']); end;
        
        Var0_tmp        = eval(['V_Img.' is{1} '_ts']);
        
        Var1_tmp        = zeros(I0,size(Var0_tmp,2));
        Var1_tmp(idx,:) = Var0_tmp;
        Y_tmp           = reshape(Var1_tmp,[X0 Y0 Z0 size(Var1_tmp,2)]);
        
        V_tmp                   = V1;
        V_tmp.hdr.dime.dim(2:5) = [X0 Y0 Z0 size(Var1_tmp,2)];
        V_tmp.img               = Y_tmp;
        save_untouch_nii(V_tmp,[DestDir is{1} '_ts.nii.gz']);
        if verbose; disp([is{1} ' saved: ' DestDir is{1} '_ts.nii.gz']); end; 
        clear *_tmp
        
        Var0_tmp      = eval(['V_Img.' is{1}]);
        Var1_tmp      = zeros(I0,size(Var0_tmp,2));
        Var1_tmp(idx) = Var0_tmp;
        Y_tmp         = flipud(reshape(Var1_tmp,[X0 Y0 Z0])); %flip back here because save_nii flips it!
        nii_tmp       = make_nii(sum(Y_tmp,4),[2,2,2],[0,0,0],64,['3D image of ' is{1}]);
        save_nii(nii_tmp,[DestDir is{1} '.nii.gz'])
        if verbose; disp([is{1} ' saved: ' DestDir is{1} '.nii.gz']); end;
        clear *_tmp
    end
    clear V_Img;
else
    if verbose
        disp('-Variance images will NOT be saved:')
        disp('-- Either destination directory was not set OR the input is not a nifti.')
    end
    clear V_Img;
end
%% Global - Res (SED vars)
V.w_Avar  = sum(V.Avar_ts);
V.w_Dvar  = sum(V.Dvar_ts);
V.w_Svar  = sum(V.Svar_ts);
V.w_Evar  = sum(V.Evar_ts);
%Global
V.g_Avar_ts = B.Ybar.^2;
V.g_Dvar_ts = B.Dbar.^2./4;
V.g_Svar_ts = B.Sbar.^2./4;
V.g_Evar_ts = B.Ybar([1,T0]).^2./2;

% Global ts (Just for vis)
V.g_Ats=B.Ybar;
V.g_Dts=B.Dbar./2;
V.g_Sts=B.Sbar./2;

V.g_Avar  = sum(V.g_Avar_ts);
V.g_Dvar  = sum(V.g_Dvar_ts);
V.g_Svar  = sum(V.g_Svar_ts);
V.g_Evar  = sum(V.g_Evar_ts);

%V.g_Avar  = sum(B.Ybar.^2);
%V.g_Dvar  = sum(B.Dbar.^2)./4;
%V.g_Svar  = sum(B.Sbar.^2)./4;
%V.g_Evar  = sum(B.Ybar([1,T0]).^2)./2;
%Non-Global

V.ng_Avar_ts = mean((Y-repmat(B.Ybar,[I1,1])).^2);
V.ng_Dvar_ts = mean((D-repmat(B.Dbar,[I1,1])).^2)./4;
V.ng_Svar_ts = mean((S-repmat(B.Sbar,[I1,1])).^2)./4;
V.ng_Evar_ts = mean([(Yhead-B.Y1bar).^2,(Ytail-B.Ytbar).^2])./2;

V.ng_Avar = sum(V.ng_Avar_ts);
V.ng_Dvar = sum(V.ng_Dvar_ts);
V.ng_Svar = sum(V.ng_Svar_ts);
V.ng_Evar = sum(V.ng_Evar_ts);

V.GrandMean_Untouched  = mean(mvY_Untouched);
V.GrandMean_NormInt    = mean(mvY_NormInt);
V.GrandMean_Demeaned   = mean(mvY_Demeaned);
V.GranMean_WholeBrain  = mean(mvY_WholeImage);
%V.ng_Avar = sum(mean((Y-repmat(B.Ybar,[I1,1])).^2));
%V.ng_Dvar = sum(mean((D-repmat(B.Dbar,[I1,1])).^2))./4;
%V.ng_Svar = sum(mean((S-repmat(B.Sbar,[I1,1])).^2))./4;
%V.ng_Evar = sum(mean([(Yhead-B.Y1bar).^2,(Ytail-B.Ytbar).^2]))./2;

% Sanity Chek - The moment of truth!
% gvars             = V.g_Dvar+V.g_Svar+V.g_Evar;
% rgvars            = V.rg_Dvar+V.rg_Svar+V.rg_Evar;
% WholeWholeVar_Test=gvars+rgvars;
% %assert(WholeWholeVar==WholeWholeVar_Test,'VarDecomp failed')
% disp(['WholeVar= ' num2str(V.w_Avar) ' , sum of decomp var= ' num2str(WholeWholeVar_Test)])
%% SED ANOVA table
SS      = I1*[V.w_Avar,V.w_Dvar,V.w_Svar,V.w_Evar,...
              V.g_Avar,V.g_Dvar,V.g_Svar,V.g_Evar];
MS      = SS/I1/T0;
RMS     = sqrt(MS);
Prntg   = RMS.^2./RMS(1).^2*100;
Expd    = [1,(T0-1)/T0/2,(T0-1)/T0/2,1/T0,...
          [1,(T0-1)/T0/2,(T0-1)/T0/2,1/T0]./I1];
RelVar  = Prntg./100./Expd;

Var_Tab = [V.w_Avar,V.w_Dvar,V.w_Svar,V.w_Evar;...
    V.g_Avar,V.g_Dvar,V.g_Svar,V.g_Evar;...
    V.ng_Avar,V.ng_Dvar,V.ng_Svar,V.ng_Evar];

DSETable = array2table([MS',RMS',Prntg',RelVar'],'VariableNames',Col_labs,'RowNames',Row_labs);

if ~isempty(DestDirTable)
   writetable(DSETable,DestDirTable) 
end

if verbose
    disp('----------------------')
    disp('Sum-of-Mean-Squared (SMS) Table')
    disp(array2table(fix(Var_Tab),'VariableNames',t3_varn,'RowNames',t3_rown))
    disp('------------')
    disp(DSETable)
    disp('----------------------')
end    

%DSE ANOVE table
Stat.Labels     = Row_labs;
Stat.SS         = SS;
Stat.MS         = MS;
Stat.RMS        = RMS;
Stat.Prntg      = Prntg;
Stat.RelVar     = RelVar;
Stat.VT         = Var_Tab;

%Config
Stat.dim        = [I1 T0]; %Inter Cranial sizes
Stat.dim0       = [I0 T0]; %the 4D image initial dimensions

%Standardised measures
Stat.DeltapDvar     = (V.Dvar_ts-median(V.Dvar_ts))./mean(V.Avar_ts)*100; % This is \Delta\%D-var i.e. How much it exceeded from it is *median* normalised by A-var. 
Stat.DeltapSvar     = (V.Svar_ts-median(V.Svar_ts))./mean(V.Avar_ts)*100; % This is \Delta\%S-var i.e. How much it exceeded from it is *median* normalised by A-var. 

Stat.pDvar          =  V.Dvar_ts./mean(V.Avar_ts)*100; % & this is \%D-var. NB! *_ts is sum across *voxels*, I in nominator and denominator cancel out. 
Stat.pSvar          =  V.Svar_ts./mean(V.Avar_ts)*100; % & this is \%S-var

%Mean -- 4 sanity checks
Stat.GranMean_WholeBrain  = mean(mvY_WholeImage);
Stat.GrandMean_Untouched  = mean(mvY_Untouched);
Stat.GrandMean_NormInt    = mean(mvY_NormInt);
Stat.GrandMean_Demeaned   = mean(mvY_Demeaned); 

function gsrY=fcn_GSR(Y)
%Global Signal Regression. From FSLnets.
%For the fMRIDiag, it needs to be transposed. 
% SA, UoW, 2017

Y=Y';
mgrot=mean(Y,2); 
gsrY=Y-(mgrot*(pinv(mgrot)*Y));
gsrY=gsrY';
