function fsl_bptf(Vin,Vout,bp,TR)
% bp should be in frq

f2tr=@(f,TR) 0.5*(1/f)/TR;
system(['/usr/local/fsl/bin/fslmaths ' Vin ' -bptf ' num2str(round(f2tr(bp(1),TR),2)) ' ' num2str(round(f2tr(bp(2),TR),2)) ' ' Vout ])