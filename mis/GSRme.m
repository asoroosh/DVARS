function gsrY=GSRme(Y,T,varargin)
% Quick Global Signal Regression
%
% NOTE: Y should be a TxI matrix
%
% SA, Ox, 2018

if nargin==3; error('something is wrong mate!'); end; 

if ~exist('T','var')
    error('Second input should be length of TS for sanity check.'); 
end; 


if size(Y,1)~=T
    Y=Y'; %TxI
end

if sum(strcmpi(varargin,'GS'))
    mY = varargin{find(strcmpi(varargin,'GS'))+1};
    if size(mY,1)~=T; error('GSRme :: Make sure that the global signal is as long as the other signals, maybe transpose?'); end
else
    mY=mean(Y,2); 
end



gsrY=Y-(mY*(pinv(mY)*Y));

if size(gsrY,1)~=size(Y,1) || size(gsrY,2)~=size(Y,2)
    error('something is wrong!')
end