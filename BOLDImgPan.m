function [Y]=BOLDImgPan(V0,varargin)

% Add toolbox to open the images-------
%addpath Nifti_Util
% Input Check-------------------------
gsrflag_lab={'NoGSR'};
gsrflag=0; verbose=1; DestDir=[]; md=[]; scl=[]; SaveMe=0; Idx=[];
if sum(strcmpi(varargin,'gsrflag'))
   gsrflag      =   varargin{find(strcmpi(varargin,'gsrflag'))+1};
end
if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
end
if sum(strcmpi(varargin,'destdir'))
   DestDir      =   varargin{find(strcmpi(varargin,'destdir'))+1};
   SaveMe       =   1;
end
if sum(strcmpi(varargin,'colrng'))
   ColRng      =   varargin{find(strcmpi(varargin,'colrng'))+1};
else
    ColRng     =   [-10 10];
end
if sum(strcmpi(varargin,'fd'))
   FDts      =   varargin{find(strcmpi(varargin,'fd'))+1};
end
if sum(strcmpi(varargin,'prefix'))
   prefix    =   varargin{find(strcmpi(varargin,'prefix'))+1};
end
if sum(strcmpi(varargin,'dvars'))
   DVARS     =   varargin{find(strcmpi(varargin,'dvars'))+1};
end
if sum(strcmpi(varargin,'idx'))
   Idx     =   varargin{find(strcmpi(varargin,'idx'))+1};
end
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end
if sum(strcmpi(varargin,'handle'))
    f_hdl          =   varargin{find(strcmpi(varargin,'handle'))+1};
else
    f_hdl=figure('position',[50,500,1600,1400]); 
    hold on; box on; 
end

if ischar(V0)
    [ffpathstr,ffname,ffext]=fileparts(V0);
    if verbose; disp(['-Path to the image is: ' ffpathstr]); end;
    
    if contains(ffname,'.dtseries') || contains(ffext,'.dtseries')
        if verbose; disp(['--File is CIFTI: ' ffname ffext]); end;
        V1 = ft_read_cifti(V0);
        V2 = V1.dtseries;
        I0 = size(V2,1); T0=size(V2,2);
        Y  = V2; clear V2 V1; 
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

%Remove voxels of zeros/NaNs-----------------
nan_idx    = find(isnan(sum(Y,2)));
zeros_idx  = find(sum(Y,2)==0);
idx        = 1:I0;
idx([nan_idx;zeros_idx]) = [];
Y([nan_idx;zeros_idx],:) = [];
I1 = size(Y,1); %update number of voxels
if verbose; disp(['-Extra-cranial areas removed: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]); end;
% Intensity Normalisation--------
IntnstyScl = @(Y,md,scl) (Y./md)*scl; 
if ~isempty(scl) && isempty(md)
    md  = median(mean(Y,2));
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Normalised by ' num2str(scl) '&' num2str(md) '.']); end;
elseif ~isempty(scl) && ~isempty(md)
    assert(md==1,'4D mean in scalling cannot be anything other than 1!')
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Scaled by ' num2str(scl) '.']); end;
elseif isempty(scl) && isempty(md)    
    disp('-No normalisation/scaling has been set!')
else
    error('Something is wrong with param re intensity normalisation')
end
%Centre the data-----------------------------
mvY    =    mean(Y,2);
dmeaner=    repmat(mvY,[1,T0]);
Y      =    Y-dmeaner; clear dmeaner
if verbose; disp('-Data centred.'); end;
%Data GSRed--------------------------------ONLY FOR TEST-----------------
%fcn_GSR = @(Y) (Y'-(mean(Y',2)*(pinv(mean(Y',2))*Y')))';
%gsrflag=1;
if gsrflag 
    gsrflag_lab={'GSR'};
    Y  =    fcn_GSR(Y);
    if verbose; disp('-Data GSRed.'); end;
end

%-----Hist
% figure; hold on; 
% histogram(Y(:),50)
%---------
nsp=14;
figure(f_hdl)
f_hdl_sp0=subplot(nsp,1,[1 10]);
hold on; box on; 
colormap(f_hdl_sp0,'gray');
imagesc(Y,ColRng)
set(f_hdl_sp0,'xticklabel',[])
axis tight
%colorbar

f_hdl_sp1=subplot(nsp,1,[11 12]);
hold on; box on;
plot(FDts,'b')
line([1 T0-1],[.2 .2],'color','r')

PatchMeUp(Idx);

ylabel('FD')
set(f_hdl_sp1,'xticklabel',[])

f_hdl_sp2=subplot(nsp,1,[13 14]);
hold on; box on;
plot(DVARS,'r')
ylabel('DVARS')

PatchMeUp(Idx);
 
if SaveMe
    if exist(DestDir,'dir')~=7; mkdir(DestDir); end;
    export_fig([DestDir '/BOLDImage_' gsrflag_lab{1} '_ColRng' num2str(ColRng(2)) '_' prefix '.png']);
end

function gsrY=fcn_GSR(Y)
%Global Signal Regression
%Inspired from FSLnets
%For the fMRIDiag, it needs to be transposed. 

Y=Y';
mgrot=mean(Y,2); 
gsrY=Y-(mgrot*(pinv(mgrot)*Y));
gsrY=gsrY';

function ph=PatchMeUp(Idx)
stpjmp=1;
yyll=ylim;
for ii=1:numel(Idx)
    xtmp=[Idx(ii)-stpjmp   Idx(ii)-stpjmp   Idx(ii)+stpjmp  Idx(ii)+stpjmp];
    ytmp=[yyll(1)               yyll(2)         yyll(2)        yyll(1)    ];
    ph(ii)=patch(xtmp,ytmp,[.5 .5 .5],'FaceAlpha',0.3,'edgecolor','none');
    clear *tmp
end
return 
