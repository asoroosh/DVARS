function fMRIDiag_plot(V,varargin)
% fMRIDiag_plot(V,varargin)
%
% Draws main var components (+ non-global comps) with BOLD intensity for
% further diagnosis. 
%
%
%   SA, 2017, UoW
%   srafyouni@gmail.com
%
%

FDflag  = 0;  AbsMovflag = 0;  BOLDFlag  = 0;  Idx     = [];
md      = []; scl        = []; verbose   = 1;  gsrflag = 0;
nsp     = 14; lw         = 2;  lfs       = 12;
TickScaler = 1; Thickness = 1; GrandMean = 100;
if sum(strcmpi(varargin,'gsrflag'))
   gsrflag      =   varargin{find(strcmpi(varargin,'gsrflag'))+1};
end
if sum(strcmpi(varargin,'GrandMean'))
   GrandMean      =   varargin{find(strcmpi(varargin,'GrandMean'))+1};
end
if sum(strcmpi(varargin,'fd'))
   FDts      =   varargin{find(strcmpi(varargin,'fd'))+1};
   FDflag    =   1; 
end
% if sum(strcmpi(varargin,'prefix'))
%    prefix    =   varargin{find(strcmpi(varargin,'prefix'))+1};
% end
if sum(strcmpi(varargin,'AbsMov'))
   AbsMov     =   varargin{find(strcmpi(varargin,'AbsMov'))+1};
   AbsMovflag =   1;
end
if sum(strcmpi(varargin,'idx'))
   Idx     =   varargin{find(strcmpi(varargin,'idx'))+1};
end
if sum(strcmpi(varargin,'BOLD'))
   Y        =   varargin{find(strcmpi(varargin,'BOLD'))+1};
   BOLDFlag = 1;
   nsp      = 20;
end
if sum(strcmpi(varargin,'handle'))
    f_hdl          =   varargin{find(strcmpi(varargin,'handle'))+1};
else
    f_hdl=figure('position',[50,500,1600,1400]); 
    hold on; box on; 
end
if sum(strcmpi(varargin,'colrng'))
   ColRng      =   varargin{find(strcmpi(varargin,'colrng'))+1};
else
    ColRng     =   [-10 10];
end
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
if sum(strcmpi(varargin,'TickScaler'))
   TickScaler          =   varargin{find(strcmpi(varargin,'TickScaler'))+1};
end
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end
if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
end

if sum(strcmpi(varargin,'linewidth'))
   lw      =   varargin{find(strcmpi(varargin,'linewidth'))+1};   
end

if sum(strcmpi(varargin,'fontsize'))
   lfs     =   varargin{find(strcmpi(varargin,'fontsize'))+1};
end
if sum(strcmpi(varargin,'Thick'))
   Thickness      =   varargin{find(strcmpi(varargin,'Thick'))+1};
end
%###################################################################################

Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
%FDcol=Col(2,:); % Red (/orange!)
FDcol=[.5 .5 .5];
Dcol=Col(1,:); % Blue
Scol=Col(3,:); % Yellow
Ecol=Col(4,:); % Purple

T=length(V.Avar_ts);
Time=1:T;
ds_bnd=1:T-1;

hTime=(1:(T-1))+0.5;
eTime=[1 T];

%###################################################################################

figure(f_hdl)

%---------------------------Allvar
% sph0=subplot(nsp,1,[1 2]);
% hold on; box on;
% plot(Time,sqrt(V.Avar_ts(Time)),'color',Acol,'linestyle','-','linewidth',lw-.5)
% ylabel('All','fontsize',lfs);
% axis tight;
% PatchMeUp(Idx);
% set(sph0,'ygrid','on','xticklabel',[])
% axis tight

%---------------------------FD%---------------------------
if  FDflag
    sph0=subplot(nsp,1,[1 2]);
    hold on; box on;
    yyaxis left
        plot(hTime,FDts,'color','k','linestyle','-','linewidth',lw-0.5)
        plot(hTime,ones(1,T-1)*0.5,'linewidth',lw,'linestyle','-.','color','r')
        plot(hTime,ones(1,T-1)*0.2,'linewidth',lw,'linestyle','-.','color','r')
        
        if max(FDts)>0.6
            ylim([0 max(FDts)+0.1]);
        else
            ylim([0 0.6]);
        end
        
        ylabel('FD (mm)','fontsize',lfs,'interpreter','latex','color','k');
        %axis tight;
        
        PatchMeUp(Idx,Thickness);
        
        set(sph0,'ygrid','on','xticklabel',[],'xlim',[1 T],'ytick',[0.2 0.5],'ycolor','k')
        if AbsMovflag
            yyaxis right
                plot(Time,AbsMov(:,1),'color',FDcol+[0.1,0.3,0.3],'linestyle','-','linewidth',lw-0.9)
                plot(Time,AbsMov(:,2),'color',FDcol+[0.1,0.5,0.5],'linestyle','-','linewidth',lw-0.9)
                axis tight
                ylabel('|D| (mm)','fontsize',lfs,'interpreter','latex');
        else
            yyaxis right
            set(sph0,'yticklabel',[],'ycolor','k')
        end
end

%---------------------------Whole%---------------------------
sph1=subplot(nsp,1,[3 11]); 
hold on; box on;
%title('DSE Variance Decomposition (RMS units)','fontsize',13)
yyaxis(sph1,'left')
    line(Time,sqrt(V.Avar_ts),'LineStyle','-','linewidth',lw,'color',Acol)
    line(Time,ones(1,T).*mean(sqrt(V.Avar_ts)),'LineStyle',':','linewidth',.5,'color',Acol)
    
    line(hTime,sqrt(V.Dvar_ts),'LineStyle','-','linewidth',lw,'color',Dcol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Dvar_ts)),'LineStyle',':','linewidth',.5,'color',Dcol)
    
    line(hTime,sqrt(V.Svar_ts),'LineStyle','-','linewidth',lw,'color',Scol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Svar_ts)),'LineStyle',':','linewidth',.5,'color',Scol)
    
    line(eTime,sqrt(V.Evar_ts),'LineStyle','none','Marker','o','markerfacecolor',Ecol,'linewidth',3,'color',Ecol)
    ylabel('$\sqrt{\mathrm{Variance}}$','fontsize',lfs,'interpreter','latex')  
    YLim2=ylim.^2/mean(V.Avar_ts)*100;
    set(sph1,'ycolor','k','xlim',[1 T])
yyaxis(sph1,'right')
    YTick2=PrettyTicks(YLim2,TickScaler); YTick=sqrt(YTick2);
    set(sph1,'XTick',[],'Ylim',sqrt(YLim2),'YTick',sqrt(YTick2),'YtickLabel',num2str([YTick2']));
    ylabel('\% of A-var','fontsize',lfs,'interpreter','latex')
    %set(sph2,'ygrid','on')
    h=abline('h',YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(sph1,'ycolor','k','xlim',[1 T])
    
PatchMeUp(Idx,Thickness);

%---------------------------Global%---------------------------
sph2=subplot(nsp,1,[12 13]); 
hold on; box on;
%yyaxis(sph2,'left')
    cntrd_g_ts=V.g_Ats+GrandMean;
    plot(Time,cntrd_g_ts,'color',Acol,'linestyle','-','linewidth',lw);
    %line(hTime,ones(1,T-1).*mean(V.g_Ats+mean(V.MeanOrig)),'LineStyle','-.','linewidth',.5,'color',Acol)
    
%%%%%%%%%%%%%%%%%%%% Un-ccomment next 4 code lines if you need to see the gDvar and gSvar. %%%%%%%%%%%%%%%%%%%% 
    %plot(hTime,V.g_Dts+mean(V.MeanOrig),'color',Dcol,'linestyle','-','linewidth',lw);
    %line(hTime,ones(1,T-1).*(mean(V.g_Dts)+mean(V.MeanOrig)),'LineStyle','-.','linewidth',.5,'color',Dcol)
    
    %plot(hTime,V.g_Sts+mean(V.MeanOrig),'color',Scol,'linestyle','-','linewidth',lw);
    %line(hTime,ones(1,T-1).*(mean(V.g_Sts)+mean(V.MeanOrig)),'LineStyle','-.','linewidth',.5,'color',Scol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mx_cntrd_g_ts=max(cntrd_g_ts); mn_cntrd_g_ts=min(cntrd_g_ts);
    stps=abs(round(diff([mx_cntrd_g_ts mn_cntrd_g_ts])./3,1));
    Ytcks=round(min(cntrd_g_ts):stps:max(cntrd_g_ts),2);
    
    ylabel('A$_{Gt}$','fontsize',lfs,'interpreter','latex')
    %axis tight
    %set(sph2,'ycolor','k')
    %Ylim=ylim; Ylim=mean(Ylim)+0.5*[-1,1]*diff(Ylim)*2; ylim(Ylim)
    %ylim_tmp=ylim; dylim_tmp=(ylim_tmp-mean(V.MeanOrig)); dylims_tmp=dylim_tmp./abs(dylim_tmp);
    %YLim22=(dylim_tmp.^2/mean(V.Avar_ts));
    %YLim22=((dylim_tmp.^2/mean(V.Avar_ts))-mean(YLim22))*100;
    set(sph2,'ygrid','on','xlim',[1 T],'ycolor','k','yTick',Ytcks)
    ytickformat('%,.2f')
    
axis tight
PatchMeUp(Idx,Thickness);
%---------------------------The big dude%---------------------------
if BOLDFlag
    if ~isnumeric(Y) && size(Y,1)<=size(Y,2); error('Unknown BOLD intensity image!'); end
    I0= size(Y,1); T0 = size(Y,2);
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
        %gsrflag_lab={'GSR'};
        Y  =    fcn_GSR(Y);
        if verbose; disp('-Data GSRed.'); end;
    end
    
    
    sph3=subplot(nsp,1,[15 20]);
    hold on; box on;
    colormap(sph3,'gray');
    imagesc(Y,ColRng)
    ylabel('Voxels','fontsize',lfs,'interpreter','latex')
    set(sph3,'xticklabel',[])
    axis tight
end

xlabel('Scans','fontsize',lfs,'interpreter','latex')
set(gcf,'Color','w');

%###################################################################################
% function T=Ticks(Ys,sph)
% % For a bunch of (Y-axis) values, find default tick locations
% % Ys - cell array of vectors to plot
% f=figure('visible','off');
% plot(Ys{1})
% hold on
% for i=2:length(Ys)
%   plot(sph,Ys{i})
% end
% hold off
% T=get(gca,'Ytick');
% close(f)
% return
%###################################################################################

function T=PrettyTicks(Lim,varargin)
% For a given axis limit, find pretty tick spacing; assumes 50 is always
% in the plot (i.e. that rounded integers are always appropriate)
% Ylim - Y axis limts
%
%   TEN & SA, 2017, UoW
%   srafyouni@gmail.com
%
MinTick=3;  % Minimum number of tick locations
if ~isempty(varargin)
    TickSp=[15 5 2.5 1 0.5 0.2]./varargin{1};
elseif isempty(varargin)
    TickSp=[15 5 2.5 1 0.5 0.2];
end
ts=0;
T=[];
while length(T)<MinTick
  ts=ts+1;
  if ts>length(TickSp)
    break
  end
  TS=TickSp(ts);
  T = ceil(Lim(1)/TS)*TS : TS : floor(Lim(2)/TS)*TS;
end

return
%###################################################################################
function h=abline(a,b,varargin)
% FORMAT h = abline(a,b,...)
% Plots y=a+b*x in dotted line
% FORMAT h = abline('h',y,...)
% Plots a horizontal line at y; y can be a vector, & then multiple lines plotted
% FORMAT h = abline('v',x,...)
% Plots a vertical line at x; x can be a vector, & then multiple lines plotted
%
% ...  Other graphics options, e.g. "'LineStyle','-'" or "'LineWidth',2" or
%      "'color',[1 0 0]",  etc
%
% Like Splus' abline.  Line is plotted and then moved behind all other
% points on the graph.
%
% $Id: abline.m,v 1.1 2013/06/04 10:38:11 nichols Exp $

if (nargin==2) && ischar(a)
  a = lower(a);
else

  if (nargin<1)
    a = 0;
  end
  if (nargin<2)
    b = 0;
  end
  
end

XX=get(gca,'Xlim');
YY=get(gca,'Ylim');
h_exist = get(gca,'children');

g = [];
if ischar(a) && (a=='h')

  for i=1:length(b)
    g=[g;line(XX,[b(i) b(i)],'LineStyle',':',varargin{:})];
  end

elseif ischar(a) && (a=='v')

  for i=1:length(b)
    g=[g;line([b(i) b(i)],YY,'LineStyle',':',varargin{:})];
  end

else

  g=line(XX,a+b*XX,'LineStyle',':',varargin{:});

end

uistack(h_exist,'top');

if (nargout>0)
  h=g;
end

set(gcf,'color','w');
return
%###################################################################################
function ph=PatchMeUp(Idx,varargin)
%   Draws a patch across the significantly identified scans on var plots
%
%   SA, 2017, UoW
%   srafyouni@gmail.com
if nargin == 1
    stpjmp=1;
elseif nargin == 2
    stpjmp=varargin{1};
end

yyll=ylim;
for ii=1:numel(Idx)
    xtmp=[Idx(ii)-stpjmp   Idx(ii)-stpjmp   Idx(ii)+stpjmp  Idx(ii)+stpjmp];
    ytmp=[yyll(1)               yyll(2)         yyll(2)        yyll(1)    ];
    ph(ii)=patch(xtmp,ytmp,[.5 .5 .5],'FaceAlpha',0.3,'edgecolor','none');
    clear *tmp
end
return 
%###################################################################################
function gsrY=fcn_GSR(Y)
%Global Signal Regression
%Inspired from FSLnets
%For the fMRIDiag, it needs to be transposed. 
%
%   SA, 2017, UoW
%   srafyouni@gmail.com
%
Y=Y';
mgrot=mean(Y,2); 
gsrY=Y-(mgrot*(pinv(mgrot)*Y));
gsrY=gsrY';