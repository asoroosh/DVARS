function [Y,Stat]=CleanNIFTI_fsl(V0,varargin)
%[Y,Stat]=CleanNIFTI_fsl(V0,varargin)
% Nomalise, demean and clean (remove zeros and NaNs) from fMRI data. 
%
% Output Y is a IxT matrix
% Trigger 'DestDir' with the destination to save the NIFTI file. 
%
% Main difference from CleanNIFTI is that this function uses FSL and is
% operational on Octave.
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________


Steps=[]; verbose=1; DestDir=[]; md=[]; scl=[];
datatype  = 'f';
voxelsize = [2 2 2 1];

if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1}; 
end
%
if sum(strcmpi(varargin,'destdir'))
   DestDir      =   varargin{find(strcmpi(varargin,'destdir'))+1};
end

if sum(strcmpi(varargin,'voxelsize'))
   voxelsize      =   varargin{find(strcmpi(varargin,'voxelsize'))+1};
end

% untested, so I leave it as it is for now
% if sum(strcmpi(varargin,'bptf'))
%     bp      =   varargin{find(strcmpi(varargin,'bptf'))+1};
%     if numel(bp<3); error('should be: [lowerbound upperbound TR]'); end;
% else
%     bp=[];
% end
%
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
%
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end
%
if sum(strcmpi(varargin,'demean'))
    dmflag          =   1;
else
    if verbose; warning('the demean flag is off!'); end;
    dmflag          =   0;
end

%sort out the FSL directories ==============================
fsldir = getenv('FSLDIR');
if isempty(fsldir) 
    error('CleanNIFTI_fsl:: I can not find FSL directory!'); 
end
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;


if ischar(V0)
    [ffpathstr,ffname,ffext]=fileparts(V0);
    if verbose; disp(['-Path to the image is: ' ffpathstr]); end;
    
    if contains(ffname,'.dtseries') || contains(ffext,'.dtseries')
        if verbose; disp(['--File is CIFTI: ' ffname ffext]); end;
        error('We can not support CIFTI files yet.')
%         V1=ft_read_cifti(V0);
%         V2=V1.dtseries;
%         I0=size(V2,1); T0=size(V2,2);
%         Y=V2; clear V2 V1; 
    elseif ~contains(ffname,'.dtseries') || contains(ffname,'.nii') 
        if verbose; disp(['--File is NIFTI: ' ffname ffext]); end;
        
        [V2,dims,voxelsize] = read_avw(V0);
        X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
        if sum(ismember(dims,[X0,Y0,Z0,T0]))~=4
            error('CleanNIFTI_fsl:: something is wrong with the dimensions!'); 
        end
        
        I0 = prod([X0,Y0,Z0]);
        Y  = reshape(V2,[I0,T0]); clear V2;
    else
        error('Unknown input image.')
    end
    
    if verbose; disp('-Image loaded.'); end;
elseif isnumeric(V0) && size(V0,1)>size(V0,2)
    if verbose; disp('-Input is a Matrix.'); end;
    Y = double(V0);  
    I0= size(Y,1); T0 = size(Y,2);
elseif isnumeric(V0) && size(V0,1)<=size(V0,2) && size(V0,3)==1
    if verbose; disp('-Input is a Matrix.'); end;
    error('Check the input, matrix should be in form of IxT, where I=XxYxZ!');    
elseif isnumeric(V0) && size(V0,3)>1 && size(V0,4)>1
    disp(['-Input is a 4D matrix.'])
    X0 = size(V0,1); Y0 = size(V0,2); Z0 = size(V0,3); T0 = size(V0,4);
    I0 = prod([X0,Y0,Z0]);
    Y  = reshape(V0,[I0,T0]); clear V2;
    V1 = [];
end

%Y = double(Y);%to work with int 16bit as well.

%------------------------------------------------------------------------
%Remove voxels of zeros/NaNs---------------------------------------------------
nan_idx    = find(isnan(sum(Y,2)));
zeros_idx  = find(sum(Y,2)==0);
idx        = 1:I0;
idx([nan_idx;zeros_idx]) = [];
Y([nan_idx;zeros_idx],:) = [];
I1 = size(Y,1); %update number of voxels

if verbose; disp(['-Extra-cranial areas removed: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]); end;
Steps=[Steps 'CLEANED_'];

Stat.GlobalMeanSignal = mean(Y);
Stat.OrigDim     = [I0 T0];
Stat.CleanedDim  = [I1 T0];
Stat.Removables  = [nan_idx;zeros_idx];
OrigMean = mean(Y,2);
Stat.OrigMean    = OrigMean;
Stat.ImgDim = [X0 Y0 Z0 T0];
Stat.voxelsize = voxelsize;
Stat.datatype  = datatype;
%------------------------------------------------------------------------
% Intensity Normalisation------------------------------------------------------
IntnstyScl = @(Y,md,scl) (Y./md).*scl; 
if ~isempty(scl) && isempty(md)
    md  = median(mean(Y,2)); %NB median of the mean image.
    %md  = mean(mean(Y,2)); %NB *mean* of the mean image.
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Normalised by ' num2str(scl) '&' num2str(md) '.']); end;
    Steps = [Steps 'NORM_'];
elseif ~isempty(scl) && ~isempty(md)
    assert(md==1,'4D mean in scalling cannot be anything other than 1!')
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Scaled by ' num2str(scl) '.']); end;
    Steps = [Steps 'SCALE_'];
elseif isempty(scl) && isempty(md)    
    if verbose; disp('-No normalisation/scaling has been set!'); end;
else
    error('Something is wrong with param re: intensity normalisation')
end
%------------------------------------------------------------------------
%Centre the data-----------------------------
if dmflag
    mvY_NormInt      =    mean(Y,2); %later will be used as grand mean! don't touch it!
    dmeaner          =    repmat(mvY_NormInt,[1,T0]);
    Y                =    Y-dmeaner; clear dmeaner
    DemeanedMen      = mean(Y,2);
    Stat.DemeanedMen = DemeanedMen;
    Steps            = [Steps 'DEMEANED_'];
end
%------------------------------------------------------------------------
%Save the image-----------------------------
if ~isempty(DestDir) && ischar(V0)
    %if ~any(strfind(path,'spm')); warning('**SPM has not been added to the path!**'); end;
    [spathstr,sname,stext]=fileparts(DestDir);
    if exist(spathstr,'dir')~=7; mkdir(spathstr); end;
                    
        Y_tmp        = zeros(I0,T0);
        Y_tmp(idx,:) = Y;
        Y_tmp        = reshape(Y_tmp,[X0 Y0 Z0 T0]);
        
        V_tmp                   = V1;
        V_tmp.hdr.dime.dim(2:5) = [X0 Y0 Z0 T0];
        save_untouch_nii(Y_tmp,[spathstr '/' sname stext],datatype,voxelsize);
        if verbose; disp(['image saved: ' spathstr '/' sname stext]); end; 
        clear *_tmp        
    clear V_Img;
else
    if verbose
        disp('-the images will NOT be saved:')
        disp('-- Either destination directory was not set OR the input is not a nifti.')
    end
    clear V_Img;
end

%------------------------------------------------------------------------ 
%Temporal filtering-----------------------------
% if ~isempty(bp)
%     disp(['bandpass filter of lowerbound ' num2str(bp(1)) ' & higherbound ' um2str(bp(2)) ' is runing!'])
%     fsl_bptf(DestDir,DestDir,[bp(1) bp(2),TR])
%     Steps=[Steps 'BPTF'];
% end
Stat.OutDir=DestDir;
Stat.Steps=Steps;

% function fsl_bptf(Vin,Vout,bp,TR)
% % bp should be in frq
% f2tr=@(f,TR) 0.5*(1/f)/TR;
% system(['/usr/local/fsl/bin/fslmaths ' Vin ' -bptf ' num2str(round(f2tr(bp(1),TR),2)) ' ' num2str(round(f2tr(bp(2),TR),2)) ' ' Vout ])