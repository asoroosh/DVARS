function gsrY=GSRme(ts)

mgrot=mean(ts,2); 
gsrY=ts-(mgrot*(pinv(mgrot)*ts));

if size(gsrY,1)~=size(ts,1) || size(gsrY,2)~=size(ts,2)
    error('something is wrong!')
end