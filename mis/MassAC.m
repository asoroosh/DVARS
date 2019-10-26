function [ACv,CLv]=MassAC(Y,nlg)
% [ACv,CLv] = MassAC(Y,nlg)
% Fast multi-lag autocorrelation estimation for large matrices. Exploit the
% matrix based features of matlab.
%
%   NB: You need MATLAB 2016b< to run the code!
%
%%%INPUTS
%   Y  :    Matrix of TxI size
%   nlg:    number of lags, default: N-1       
%
%%%OUTPUTS
%   ACv:    Matrix of AC coefficients #lags x #Voxels
%   CLv:    Correlation length of size 1 x #Voxels
%   
%   SA - 2017
%

if size(Y,1)>size(Y,2)
    error('Y should be TxI')
end

if ~exist('nlg','var')
    nlg=size(Y,1)-1;
end

I0      = size(Y,2);
ndpr    = size(Y,1);

%----------------------------------------------------

Ybar=mean(Y);
c0=(1/(ndpr-1))*sum((Y-Ybar).^2);

for lg=1:nlg
    if ~mod(lg,10); disp(['on lag: ' num2str(lg)]); end;
    
    Y0=Y(1:end-lg,:);
    Y1=Y(lg+1:end,:);
    
    ck=(1/(ndpr-1))*sum((Y0-Ybar).*(Y1-Ybar));
    
    ACv(lg,:)=ck./c0;
    
    clear ck Y0 Y1 stdY0 stdY1
end
if nlg==1
    CLv=ACv; %because it is how it is :)
else
    CLv=1+sum(ACv.^2); %because we didn't take into acount the lag-0!
end
