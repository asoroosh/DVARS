function [idx,val,Dv]=DVARS_IQR(YDv,varargin)
%
% [idx,val]=DVARS_IQR(YDv,varargin)
%
%   Filter (or Calculate and then filter) the DVARS as proposed in fsl_motion_outliers. 
%
%   For filtering the 
%   [idx,val]=DVARS_IQR(DVARS,'DVARS')
%
%   For calculating DVARS and then Filter it. 
%   [idx,val]=DVARS_IQR(Y,'BOLD')
%____________________________________________
%
%   SA, NISOx, 2017
%   srafyouni@gmail.com

if sum(strcmpi(varargin,'DVARS'))
    assert(ismember(1,size(YDv)),'Input is not DVARS!')
    Dv=YDv;
elseif sum(strcmpi(varargin,'BOLD'))
    assert(~ismember(1,size(YDv)) && size(YDv,1)>size(YDv,2),'Input is not BOLD or the input should be transposed!')
    Y=YDv;
    I=size(Y,1);
    DY    = diff(Y,1,2);
    Dv = sqrt(sum(DY.^2)./I); clear DY Y; %save me some memory!!
end
%
Q1 = quantile(Dv,0.25);
Q3 = quantile(Dv,0.75);
IQR = 1.5*(Q3-Q1);
UpperBound = Q3 + IQR;
LowerBound = Q1 - IQR;
idx=find( Dv>UpperBound | Dv<LowerBound);
val=Dv(idx);