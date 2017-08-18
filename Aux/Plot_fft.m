function [Y_ft,freq,Pxx]=Plot_fft(Y_t,L,TR,varargin)
% SA, NISOx, 2017

lfts=12; %font size for labels

if sum(strcmpi(varargin,'subplot'))
   hndl  = varargin{find(strcmpi(varargin,'subplot'))+1};
   subplot(hndl)
elseif sum(strcmpi(varargin,'figure'))
   hndl  = varargin{find(strcmpi(varargin,'figure'))+1};
   figure(hndl)
else
    hndl = figure; hold on; box on;
end

if size(Y_t,2)~=L
    Y_t=Y_t';
end

Fs      = 1/TR;              %sampling frequency

nfft    = 2^nextpow2(2*L-1);
Y_ft     = fft(Y_t,nfft,2); % not now, but later, fix the frequency resolution/length of the fft

Pxx=1/(L*Fs)*abs(Y_ft(:,1:L/2+1)).^2;
Pxx(:,2:end-1) = 2*Pxx(:,2:end-1);

freq = 0:Fs/L:Fs/2;

% figure; hold on; 
% subplot(2,1,1)
% hold on; box on; grid on;
% plot(freq,Pxx,'color',[.5 .5 .5 .5]) 
% plot(freq,mean(Pxx),'r-.','linewidth',1.4)
% xlabel('Hz')
% ylabel('|P1(f)|')

%subplot(2,1,2)
hold on; box on; grid on;
plot(freq,10*log10(Pxx),'color',[.5 .5 .5 .2]) 
plot(freq,mean(10*log10(Pxx)),'r-','linewidth',1.4)
xlabel('Hz','fontsize',lfts,'interpreter','latex'); 
ylabel('dB/Hz','fontsize',lfts,'interpreter','latex');
%-------------------------------


