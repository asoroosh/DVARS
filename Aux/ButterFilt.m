function Y=ButterFilt(Y,Ord,bnd_freq,TR)
% Butterworth filtering of the multi-dim data
% Input Y should be in IxT form. 
%
% SA, NISOx.org, 2017
%

% if size(Y,1)<size(Y,2)
%     error('The diminsions should be as IxT.')
% end

if isempty(Ord); Ord=4;         end; %just to be safe!
if isempty(bnd_freq); bnd_freq=[.01 .1];  end;

srate = 1/TR;
bnd_radsmpl   = bnd_freq./(srate/2);

if  bnd_radsmpl(1)<0 && numel(bnd_radsmpl)==2
    disp(['-Highpass filter (Butterworth) of order ' num2str(Ord) ' & cutoff: ' num2str(bnd_freq(2)) ' Hz'])
    [b,a] = butter(Ord,bnd_radsmpl(2),'high');
    Y     = filter(b,a,Y,[],2);
elseif bnd_radsmpl(2)<0 && numel(bnd_radsmpl)==2
    disp(['-Lowepass filter (Butterworth) of order ' num2str(Ord) ' & cutoff: ' num2str(bnd_freq(1)) ' Hz'])
    [b,a] = butter(Ord,bnd_radsmpl(1),'low');
    Y     = filter(b,a,Y,[],2);    
elseif and(bnd_radsmpl(1)>0,bnd_radsmpl(2)>0) && diff(bnd_radsmpl)>0 && numel(bnd_radsmpl)==2
    disp(['-Bandpass filter (Butterworth) of order ' num2str(Ord) ' & band : [' num2str(bnd_freq(1)) ' ' num2str(bnd_freq(2)) '] Hz'])
    [b,a] = butter(Ord,bnd_radsmpl);
    Y     = filter(b,a,Y,[],2);     
else
    error('The band is missing or maldefined!')
end