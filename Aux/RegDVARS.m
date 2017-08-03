function dv_Y=RegDVARS(Y,rgsrs,varargin)

md          = [];   scl         = [];

% if size(Y,1)<size(Y,2)
%     Y=Y';
%     warning(['dude reshaped to IxT form: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]);
% end

%
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end
%
if sum(strcmpi(varargin,'verbose'))
    verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
else
    verbose     = 1;
end
%--------------------------------------------------------------------------
IntnstyScl = @(Y,md,scl) (Y./md).*scl; 

if ~isempty(scl) && isempty(md)
    md  = median(mean(Y,2)); %NB median of the mean image
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
mvY     =    mean(Y,2);
%size(mvY)
dmeaner =    repmat(mvY,[1,size(Y,2)]);
Y       =    Y-dmeaner; clear dmeaner
%--------------------------------------------------------------------------
dv_Y=Y'-(rgsrs*(pinv(rgsrs)*Y'));