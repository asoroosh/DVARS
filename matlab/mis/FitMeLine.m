function [xl,yl]=FitMeLine(X,Y,Ord)
% Quickly but nicely fit a line of chosen order!
% SA-2017
for o=1:numel(Ord)
    xl_tmp=linspace(min(X),max(X),200);
    xl(o,:)=xl_tmp;
    cf_tmp=polyfit(X, Y, Ord(o));
    yl(o,:)=polyval(cf_tmp,xl_tmp);
    clear *_tmp
end