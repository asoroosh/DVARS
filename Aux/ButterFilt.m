function Y=ButterFilt(Y,Ord,bnd,TR)
% Butterworth filtering of the multi-dim data
%
% SA, NISOx.org, 2017
%
if size(Y,1)<size(Y,2)
    error('The diminsions should be as IxT.')
end

if isempty(Ord); Ord=4;         end;
if isempty(bnd); bnd=[.01 .1];  end;

srate = 1/TR;
bnd   = bnd./(srate/2)

if  bnd(1)<0 && numel(bnd)==2
    disp(['-Highpass filter.'])
    [b,a] = butter(Ord,bnd(2),'high');
    Y     = filter(b,a,Y,[],2);
elseif bnd(2)<0 && numel(bnd)==2
    disp(['-Lowpass filter.'])
    [b,a] = butter(Ord,bnd(1),'low');
    Y     = filter(b,a,Y,[],2);    
elseif and(bnd(1)>0,bnd(2)>0) && diff(bnd)>0 && numel(bnd)==2
    disp(['-Bandpass filter.'])
    [b,a] = butter(Ord,bnd);
    Y     = filter(b,a,Y,[],2);     
else
    error('The band is missing or maldefined!')
end