function gsrY=GSRme(Y,T)
% Quick Global Signal Regression
%
% NOTE: Y should be a TxI matrix
%
% SA, Ox, 2018

if size(Y,1)~=T
    Y=Y'; %TxI
end

mY=mean(Y,2); 
gsrY=Y-(mY*(pinv(mY)*Y));

if size(gsrY,1)~=size(Y,1) || size(gsrY,2)~=size(Y,2)
    error('something is wrong!')
end