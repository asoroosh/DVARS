function [DVARS,Stat]=DVARSCalc(V0,varargin)
%[DVARS,Stat]=DVARSCalc(V0)
% Statistical inference on DVars component to identify corrupted scans. 
%
%%%%INPUTS:
%
%   V0:             Can be both (1) a string indicating the path to the 
%                   nifti file OR a numerical matrix of size IxT. Where I 
%                   is number of voxels (I=Nx x Ny x Nz) and T is number of
%                   data-points.
%
%   Following parameters are optional:
%
%   'TestMethod':   Should be followed by 'Z' for Z-test and 'X2' for 
%                   Chi^2 test [default:'Z'].
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'TestMethod','X2')
%
%   'VarType':      Method for robust estimate of variance. It can be either 
%                   'IQR' for full IQR or 'hIQR' for half-IQR. 
%                   [default:'hIQR']
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'VarType','IQR')
%
%   'MeanType':     Method for robust estimate of expected value. The value
%                   should be a digit corresponding to the order of 
%                   following methods [default:'median']: 
%                   'sig2bar','sig2median','median','sigbar2','xbar'.
%                   For example: MeanType=3 means the method to estimate
%                   robust expected value is empirical median. 
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'MeanType',3)
%
%   'TransPower':   Power of transformation [default:1/3]
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'TestMethod','X2','TransPower',1/3)
%
%   'RDVARS':       By passing this arg, the function generates the
%                   relative DVARS (RDVARS). NB! this might take a while
%                   due to robust estimate of autocorrelation. 
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'RDVARS')
%
%   'verbose':      Set to 1 if you need the log of runing code [default:1]
%                   e.g. [DVARS,Stat]=DVARSCalc(V0,'verbose',1)
%
%   'Norm'          Intensity normalisation to a given scale.
%                   e.g.: [V,Stat]=DSEvars(V0,'Norm',100)
%
%   'Scale'         Scale the intensity between the data-sets.
%                   e.g.: [V,Stat]=DSEvars(V0,'Scale',1/10)
%
%%%%OUTPUTS:
%
%   DVARS:  a vector of size 1xT-1 of classic DVARS measure
%   Stat:   a structure contains all the details of the statistical inference
%           including the standardised DVARS, pvals and further summary stats. 
%
%%%%NOTES:
%   1) It is recommended to only use time series of intra-cranial voxels as
%   inclusding the extra-cranial may inflate the variance. You can use
%   'bet' in FSL package to remove the extra-cranial areas. The scripts
%   automatically remove the zero/NaN voxels.
%   
%
%
%%%%EXAMPLE:
%
%   For iid case:
%
%   I=4e4; T=1200; Y=randn(I,T);
%   [DVARS,Stat]=DVARSCalc(Y,'VarType','hIQR','TestMethod','X2','TransPower',1/3);
%   find(Stat.pvals<0.05./(T-1)) %print corrupted DVARS data-points
%
%   For the case with simulated ouliers:
%
%   I=4e4; T=1200;
%   Y=randn(I,T);
%   Idx_OL=randi(T);
%   Y(:,Idx_OL)=Y(:,Idx_OL)+1;
%   [DVARS,Stat]=DVARSCalc(Y,'VarType','hIQR','TestMethod','X2','TransPower',1/3);
%   find(Stat.pvals<0.05./(T-1)) %print corrupted DVARS data-points
%
%%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021
%
%   Soroosh Afyouni & Thomas Nichols, UoW, Feb 2017
% 
%   https://github.com/asoroosh/DVARS
%   http://warwick.ac.uk/tenichols
%
%   Please report bugs to srafyouni@gmail.com
%

%ParCheck------------------------------------------------------------------
NDVARS_X2 = 'N/A'; NDVARS_Z = 'N/A'; RelDVARS = 'N/A';

testmeth    = 'Z';  nflag       = 0;
dd          = 1;    verbose     = 1;
WhichExpVal = 3;    WhichVar    = 3;
ACf_idx     = [];   gsrflag     = 0; 
md          = [];   scl         = [];

% Add toolbox to open the images-------------------------------------------
addpath Nifti_Util
% Input Check--------------------------------------------------------------
if sum(strcmpi(varargin,'gsrflag'))
   gsrflag      =   varargin{find(strcmpi(varargin,'gsrflag'))+1};
end
if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
end
if sum(strcmpi(varargin,'TestMethod'))
   testmeth     =   varargin{find(strcmpi(varargin,'TestMethod'))+1};
   if strcmpi(testmeth,'Z')
       WhichExpVal = 3;
   end
end
if sum(strcmpi(varargin,'transpower'))
   dd           =   varargin{find(strcmpi(varargin,'transpower'))+1};
end
if sum(strcmpi(varargin,'RDVARS'))
   %nflag        =   varargin{find(strcmpi(varargin,'RDVARS'))+1};
   nflag        =   1;
end
if sum(strcmpi(varargin,'MeanType'))
   WhichExpVal  =   varargin{find(strcmpi(varargin,'MeanType'))+1};
end
if sum(strcmpi(varargin,'VarType'))
   switch varargin{find(strcmpi(varargin,'VarType'))+1}
    case 'IQR'
        WhichVar = 2;
    case 'hIQR'
        WhichVar = 3;
    otherwise
        error('Unknown VarRobIQR! Choose between IQR and hIQR! NB case sensitive')
    end
end
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end
%--------------------------------------------------------------------------
if ischar(V0)
    [ffpathstr,ffname,ffext]=fileparts(V0);
    if verbose; disp(['-Path to the image is: ' ffpathstr]); end;
    
    if contains(ffname,'.dtseries') || contains(ffext,'.dtseries')
        if verbose; disp(['--File is CIFTI: ' ffname ffext]); end;
        V1=ft_read_cifti(V0);
        V2=V1.dtseries;
        I0=size(V2,1); T0=size(V2,2);
        Y=V2; clear V2 V1; 
    elseif ~contains(ffname,'.dtseries') || contains(ffname,'.nii') 
        if verbose; disp(['--File is NIFTI: ' ffname ffext]); end;
        V1 = load_untouch_nii(V0);
        V2 = V1.img; 
        X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
        I0 = prod([X0,Y0,Z0]);
        Y  = reshape(V2,[I0,T0]); clear V2 V1;
    else
        error('Unknown input image.')
    end
    
    if verbose; disp('-Image loaded.'); end;
elseif isnumeric(V0) && size(V0,1)>size(V0,2)
    if verbose; disp('-Input is a Matrix.'); end;
    Y = V0;
    I0= size(Y,1); T0 = size(Y,2);
elseif isnumeric(V0) && size(V0,1)<=size(V0,2)
    if verbose; disp('-Input is a Matrix.'); end;
    error('Check the input, matrix should be in form of IxT, where I=XxYxZ!');    
end

%Remove voxels of zeros/NaNs----------------------------------------------
nan_idx    = find(isnan(sum(Y,2)));
zeros_idx  = find(sum(Y,2)==0);
idx        = 1:I0;
idx([nan_idx;zeros_idx]) = [];
Y([nan_idx;zeros_idx],:) = [];
I1   = size(Y,1); %update number of voxels
if verbose; disp(['-Extra-cranial areas removed: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]); end;

mvY0 = mean(Y,2); % untouched grand mean 

% Intensity Normalisation----------------------------------------------
IntnstyScl = @(Y,md,scl) (Y./md).*scl; 

if ~isempty(scl) && isempty(md)
    md  = median(median(Y,2));
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Normalised by, scale: ' num2str(scl) ' & median: ' num2str(round(md,2)) '.']); end;
    
elseif ~isempty(scl) && ~isempty(md)
    assert(md==1,'4D mean in scalling cannot be anything other than 1!')
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Scaled by ' num2str(scl) '.']); end;
    
elseif isempty(scl) && isempty(md)    
    if verbose;  disp('-No normalisation/scaling has been set!'); end;
    
else
    error('IntnstyScl :: Something is wrong with param re intensity normalisation')
end

%Centre the data-----------------------------------------------------------
mvY    =    mean(Y,2);
dmeaner=    repmat(mvY,[1,T0]);
Y      =    Y-dmeaner; clear dmeaner
if verbose; disp(['-Data centred. Untouched Grand Mean: ' num2str(mean(mvY0)) ', Post-norm Grand Mean: ' num2str(mean(mvY))]); end;
%Data GSRed--------------------------------ONLY FOR TEST-------------------
%fcn_GSR = @(Y) Y'-(mean(Y,2)*(pinv(mean(Y,2))*Y'));
%gsrflag=1;
if gsrflag 
    Y  =    fcn_GSR(Y);
    if verbose; disp('-Data GSRed.'); end;
end
%----------------------------------------^^ONLY FOR TEST-------------------

%funcs-----
IQRsd   =   @(x) (quantile(x,0.75)-quantile(x,0.25))/1.349;
H_IQRsd =   @(x) (quantile(x,0.5)-quantile(x,0.25))/1.349*2;
%--
Zstat   =   @(x,m,s) (x-m)/s;
Zp      =   @(x,m,s) 1-normcdf(Zstat(x,m,s));
%--
X2stat  =   @(x,m,s) 2*m/s^2*x;
X2df    =   @(m,s)   2*m^2/s^2;
X2p     =   @(x,m,s) chi2cdf(X2stat(x,m,s),X2df(m,s),'upper');
%X2p     =   @(x,m,s) 1-chi2cdf(X2stat(x,m,s),X2df(m,s));

X2p0    =   @(x,m,s) (X2stat(x,m,s)-X2df(m,s))/sqrt(2*X2df(m,s));
%Relative DVARS--------------------------------------------------------
DY    = diff(Y,1,2);
DVARS = sqrt(sum(DY.^2)./I1);
if nflag
    if verbose; disp(['-Robust estimate of autocorrelation...']); end; 
    Rob_S = IQRsd(Y');
    AC    = zeros(1,I1);
    for iv=1:I1
        if (~mod(iv,10e2) && verbose); disp(['--voxel: ' num2str(iv)]); end;
        AC(iv) = madicc(Y(iv,1:end-1),Y(iv,2:end));
    end
    
    ACf_idx = isnan(AC);
    
    if any(ACf_idx)  
        AC(ACf_idx)=[]; Rob_S(ACf_idx)=[]; 
        if verbose; disp(['--AC robust estimate was failed on ' num2str(sum(ACf_idx)) ' voxels.']); end; 
    end
    
    RelDVARS     =  DVARS./(sqrt((sum(2*(1-AC).*(Rob_S.^2)))./I1));
end
%Inference-----------------------------------------------------------------
DVARS2  = mean(DY.^2);
Rob_S_D = IQRsd(DY')';

MeanNms = {'sig2bar','sig2median','median','sigbar2','xbar'};
Mn = [mean(Rob_S_D.^2),median(Rob_S_D.^2),median(DVARS2),...
    mean(Rob_S_D).^2,mean(DVARS2)];

Z   = DVARS2.^dd;
M_Z = median(Z);

VarNms = {'S2' , 'IQRd' , 'hIQRd'};
Va     = [round(var(DVARS2),2) , (1/dd*M_Z^(1/dd-1)*IQRsd(Z))^2 , (1/dd*M_Z^(1/dd-1)*H_IQRsd(Z))^2];

if verbose
    disp('-Settings: ')
    disp(['--Test Method:          ' testmeth]);
    disp(['--ExpVal method:        ' MeanNms{WhichExpVal}]); 
    disp(['--VarEst method:        ' VarNms{WhichVar}]);
    disp(['--Power Transformation: ' num2str(round(dd,2))]);
end;

switch testmeth
    case 'Z'
        M_DV2 = Mn(WhichExpVal); 
        S_DV2 = sqrt(Va(WhichVar));      
        Zval  = Zstat(DVARS2,M_DV2,S_DV2);
        Pval  = Zp(DVARS2,M_DV2,S_DV2);
        %Pval(Pval==0) = 10e-15; %There is no p-value=0!!
        nu = []; c = []; NDVARS_X20=[];
        NDVARS_Z=Zval;
    case 'X2'
        M_DV2         = Mn(WhichExpVal);             
        S_DV2         = sqrt(Va(WhichVar));       
        Pval          = X2p(DVARS2,M_DV2,S_DV2);
        Zval          = Zstat(DVARS2,M_DV2,S_DV2);
        c             = X2stat(DVARS2,M_DV2,S_DV2);
        nu            = X2df(M_DV2,S_DV2);          %Spatial EDF 
        NDVARS_X2     = -norminv(Pval);
        NDVARS_X20    = X2p0(DVARS2,M_DV2,S_DV2);
        
        %only substitute the infs, the rest is better done in matlab.
        NDVARS_X2(isinf(NDVARS_X2)) = NDVARS_X20(isinf(NDVARS_X2));
        
        
    otherwise
        error('Unknown test method!')
end

if verbose
    fprintf('\nSettings: TestMethod=%s  I=%d  T=%d \n',testmeth,I1,T0)
    disp('----Expected Values----------------------------------')
    disp(array2table(Mn,'VariableNames',MeanNms));...
    disp('----Variances----------------------------------------')
    disp(array2table(Va,'VariableNames',VarNms));...
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%OLD CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(DVARS2Smethod,'TR1')
%     S_dvars2=diff(quantile(DVARS2,[0.25 0.75]))./1.349;
%     nu=(2*E_dvars2.^2)./(S_dvars2.^2);
%     c=(2*E_dvars2)./(S_dvars2.^2);
%     critval=chi2inv(1-alp/T,nu);
%     pvals=1-chi2cdf(DVARS2*c,nu);
%     pvals_adj=pvals.*T;
%     %disp('Corrupted volumes:')
%     RemoveMe=find(DVARS2*c>=critval);
% elseif strcmp(DVARS2Smethod,'TR2')
%     YD=diff(Y');
%     S_dvars2=2*(diff(quantile(sum(YD(3:T-1,:).*YD(1:T-3,:),2)/I,[0.25 0.75]))./1.349);
%     nu=(2*E_dvars2.^2)./(S_dvars2.^2);
%     c=(2*E_dvars2)./(S_dvars2.^2);
%     critval=chi2inv(1-alp/T,nu);
%     pvals=1-chi2cdf(DVARS2*c,nu);
%     pvals_adj=pvals.*T;
%     RemoveMe=find(DVARS2*c>=critval);
% elseif strcmp(DVARS2Smethod,'TransformationChi2')
%     if d_tran==1
%         S_dvars2=diff(quantile(DVARS2,[0.25 0.75]))./1.349;
%         
%     else				       
%         Z = DVARS2.^d_tran;
%         S_Z=diff(quantile(Z,[0.25 0.75]))./1.349; %robust std of Z
%         E_dvars2=E_dvars2^(1/d_tran);
%         S_dvars2=(1/d_tran*E_dvars2^(1/d_tran-1)*S_Z);
%     end
%     c=2*E_dvars2/S_dvars2^2;
%     nu=2*E_dvars2^2/S_dvars2^2;
% 
%     %pvals=1-normcdf((DVARS2-E_dvars2)/S_dvars2);
%     pvals=1-chi2cdf(DVARS2*c,nu);
%     pvals_adj=pvals.*T;
%     RemoveMe=find(pvals_adj<0.05);
%     critval=[];
%     
% elseif strcmp(DVARS2Smethod,'TransformationNormal')
%     Z=(DVARS2./E_dvars2).^d_tran;
%     S_Z=diff(quantile(Z,[0.25 0.75]))./1.349; %robust std of Z
%     S_dvars2=sqrt((1./(E_dvars2.^2))*(((1./d_tran)-1).^2)*(S_Z.^2));
%     
%     Zs=(DVARS2.^d_tran-E_dvars2.^d_tran)./S_dvars2;
%     pvals=2*(normcdf(-abs(Zs),0,1));
%     pvals_adj=pvals.*T;
%     RemoveMe=find(pvals_adj<0.05);
%     c=[]; nu=[]; critval=[];
% else
%     disp('Robust STD method for DVARS2 has not been identified correctly!')
% end
    
%MnTrue=2*mean(SD.^2);
%VaTrue=8*mean(SD.^4)/I;

%Stat.iid.Mn=MnTrue;
%Stat.iid.Vn=VaTrue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stat.DVARS2     = DVARS2;
%Test stats
Stat.E          = Mn;
Stat.S          = sqrt(Va);
Stat.nu         = nu;
Stat.c          = c;
Stat.pvals      = Pval;
Stat.Zval       = Zval;
%Standardised 
Stat.RDVARS     = RelDVARS;
Stat.SDVARS_X2  = NDVARS_X2;  %this dude is suggested to be used!
Stat.SDVARS_X20 = NDVARS_X20;
Stat.SDVARS_Z   = NDVARS_Z;
%General info
Stat.Mu0        = mean(IQRsd(Y));
Stat.Y_MS       = mean(mean(Y.^2)); 
Stat.dim        = [I1 T0];
Stat.dim0       = [I0 T0];
Stat.RobAC_Fail = find(ACf_idx);
Stat.GrandMean  = mean(mvY);
Stat.GrandMean0 = mean(mvY0);
%Stat.NDVARS_Mu0 = (DVARS-Stat.Mu0)/sqrt(Stat.Y_MS)*100; 


Stat.Config.TestMeth    =   testmeth;
Stat.Config.PowerTrans  =   dd;
Stat.Config.WhichExpVal =   MeanNms{WhichExpVal};
Stat.Config.WhichVar    =   VarNms{WhichVar};
Stat.Config.gsrflag     =   gsrflag;
Stat.Config.Median2Norm =   md;
Stat.Config.Scale2Norm  =   scl;

function rmad = madicc(x,y)
% Median Absolute Deviation Intraclass Correlation Coefficient
%
% Impliments Median Absolute Deviation Correlation Coefficient, (as
% described in Shevlyakov & Smirnov (2011)), modified to be the
% intraclass version of correlation.  The (non-intrasclass) 
% estimate is
%     r = ( mad(Sp)^2 - mad(Sm)^2 ) / ( mad(Sp)^2 + mad(Sm)^2 )
% where
%     Sp = (x-m(x))/mad(x) + (y-m(y))/mad(y);
%     Sm = (x-m(x))/mad(x) - (y-m(y))/mad(y);
% and m() is median and mad() is the median absolute deviation,
%     mad(x) = m(abs(x-m(x)))
%
% For intraclass correlation we assume mad(x)=mad(y) and so the divisors
% cancel; further, we find a common estimate of the median mm=m([x,y])
% can compute Sp & Sm as:
%     Sp = (x-mm) + (y-mm);
%     Sm = (x-mm) - (y-mm);
%
%
% REFERENCES
%
% Shevlyakov, G., & Smirnov, P. (2011). Robust estimation of a
% correlation coefficient: an attempt of survey. Australian & New
% Zealand Journal of Statistics, 40(1), 147–156.
% 
% Kharin, Y. S., & Voloshko, V. A. (2011). Robust estimation of AR 
% coefficients under simultaneously influencing outliers and missing 
% values. Journal of Statistical Planning and Inference, 141(9),
% 3276–3288.
%
% 2014-07-08
% Thomas Nichols http://warwick.ac.uk/tenichols

I=find(all(~isnan([x(:) y(:)]),2));
if isempty(I)
  rmad=NaN;
else
  mx    = median(x(I));
  my    = median(y(I));
  Sp    = (x(I)-mx) + (y(I)-my);
  Sm    = (x(I)-mx) - (y(I)-my);
  madSp = median(abs(Sp-median(Sp)));
  madSm = median(abs(Sm-median(Sm)));
  if madSp==0 && madSm==0
    rmad = NaN;
  else
    rmad = (madSp^2 - madSm^2)/(madSp^2 + madSm^2);
  end
end


function gsrY=fcn_GSR(Y)
%Global Signal Regression
%Inspired from FSLnets
%For the fMRIDiag, it needs to be transposed. 

Y=Y';
mgrot=mean(Y,2); 
gsrY=Y-(mgrot*(pinv(mgrot)*Y));
gsrY=gsrY';
