#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:56:21 2018

@author: sorooshafyouni
"""
##########################################################
##########################################################
#DO NOT USE THIS FUNCTION!! IT IS NOT FINISHED YET! 
#USE MATLAB VERSIONS INSTEAD!
##########################################################
##########################################################

#DSE('/Users/sorooshafyouni/Desktop/HCP/135932/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz')
#print "DONE!"

#def DSE(path2img,**kwargs):
from __future__ import division
import numpy as np 
import nibabel as nib #pip install nibabel
import os.path as op
import os
import matplotlib.pyplot as plt
import scipy as sp
import scipy.io as matfile
#import json

def SaveMe2txt(Dict2Save,Where2Save):
    if not os.path.exists(Where2Save):
        os.makedirs(Where2Save)
    for varname in Dict2Save.keys():       
        np.savetxt(Where2Save + '/DSE_' + str(varname) + '.txt',Dict2Save[varname],fmt='%10.5f')

def SaveMe2Nifti(Dict2Save,Where2Save,imobj,rmvIdx):
    if not os.path.exists(Where2Save):
        os.makedirs(Where2Save)

        
    imghdr = imobj.header.copy()
    [X,Y,Z,T] = imobj.shape[0:4]
    I = np.prod((X,Y,Z))
    idx_tmp  = range(I+1)[1:I+1]
    idx_tmp0 = np.delete(idx_tmp,rmvIdx,axis=0)        
    
    for varname in Dict2Save.keys():    
        print 'writing *' + varname + '* images...'
        
        if len(np.shape(Dict2Save[varname]))==1:            
            G = Dict2Save[varname];
            Y2I = np.zeros(shape=I)
            Y2I[idx_tmp0] = G;
            I3Y = np.reshape(Y2I,(X,Y,Z))
            
        elif len(np.shape(Dict2Save[varname]))==2:            
            G = Dict2Save[varname];
            Y2I = np.zeros(shape=(I,np.shape(G)[1]))
            Y2I[idx_tmp0,0:T] = G;
            I3Y = np.reshape(Y2I,(X,Y,Z,np.shape(G)[1]))
            
        else: raise ValueError ,"Something is wrong with image dimensions."
        
        imgobj_vars = nib.nifti1.Nifti1Image(I3Y, None, header=imghdr)
        nib.save(imgobj_vars, Where2Save+'/DSE_'+ str(varname) +'.nii.gz')
        
        del(Y2I,I3Y,G) 
        
scl = 100
 #   scl = kwargs.get('Scale', None)
#d = kwargs.get('d', None)
#def DSE():

################################Read and Clean
Path2Img = '/Users/sorooshafyouni/Desktop/HCP/135932/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz'
pth, fn = os.path.split(Path2Img)

NewDir2DSE = pth + '/DSE/'

filename =  op.expandvars(Path2Img)
imobj = nib.load(filename, mmap=False)
imdat = imobj.get_data().astype(float)
####################################################################
T = imdat.shape[3]
I = np.prod(imdat.shape[0:3])
Y = np.reshape(imdat,(I,T))

#IQRsd   = lambda x: (np.percentile(x,75)-np.percentile(x,25))/1.349
#H_IQRsd = lambda x: (np.percentile(x,50)-np.percentile(x,25))/(2*1.349)
##############################
#Testings    
#Y = np.random.randn(5,5,5,550);
#mat = matfile.loadmat('TestTS.mat')
#Y = mat['Y']
#T = np.shape(Y)[3]
#I = np.prod(np.shape(Y)[0:3])
#Y = np.reshape(Y,(I,T))
#############################

print "The image has total of: " + str(I) + " voxels & " + str(T) + " time-points."
#find zeros or nans
zrIdx = np.where(~Y.any(axis=1))[0];
nanIdx = np.where(np.isnan(Y))[0]; 
rmvIdx = np.concatenate((zrIdx,nanIdx));

Y = np.delete(Y,rmvIdx,axis=0)
I1 = Y.shape[0]

print "After cleaning; " + str(I1) + " voxels & " + str(T) + " time-points."

#Intensity Normalisation 
print "Intensity Normalisation to " + str(scl)

md = np.median(np.mean(Y,axis = 1))
Y = Y/md * scl;

#Demean each time series
print "Demean along T"
mY2 = np.mean(Y,axis=1)
Y = Y-np.transpose(np.tile(mY2,(T,1)));

###################################################################
bar_Y  = np.sum(Y,axis=0)/I1; #global signal is here!
#Dvar
Dvar = Y[:,0:T-1]-Y[:,1:T] #OMG, python!
bar_Dvar  = np.sum(Dvar,axis=0)/I1;

#Svar
Svar = Y[:,0:T-1]+Y[:,1:T]
bar_Svar = np.sum(Svar,axis=0)/I1;

#Evar
Ytail = Y[:,-1]; 
Yhead = Y[:,0];
bar_Ytbar = np.sum(Ytail)/I1; 
bar_Y1bar = np.sum(Yhead)/I1;

#SED Var Images######################################################
print "Variance 4D Images"
VImg_Avar_ts = Y**2;
VImg_Dvar_ts = Dvar**2/4;
VImg_Svar_ts = Svar**2/4;
VImg_Evar_ts = np.transpose(np.array((Yhead,Ytail))**2/2);

#Averaged across time#######################################
print "Variance 3D Images"
VImg_Avar = np.mean(Y**2,axis=1);
VImg_Dvar = np.mean(Dvar**2,axis=1)/2;
VImg_Svar = np.mean(Svar**2,axis=1)/2;
VImg_Evar = np.mean(np.array((Yhead**2,Ytail**2)),axis=0); # <<<< should be checked

vimg = {'Avar_Img3':VImg_Avar,'Dvar_Img3':VImg_Dvar,
        'Svar_Img3':VImg_Svar,'Evar_Img3':VImg_Evar,
        'Avar_Img4':VImg_Avar_ts,'Dvar_Img4':VImg_Dvar_ts,
        'Svar_Img4':VImg_Svar_ts,'Evar_Img4':VImg_Evar_ts};

SaveMe2Nifti(vimg,NewDir2DSE+"/VarImg",imobj,rmvIdx)
         
#DSE Time series -- averaged across I#######################
print "DSE time series"
V_Avar_ts = np.mean(VImg_Avar_ts,axis=0);
V_Dvar_ts = np.mean(VImg_Dvar_ts,axis=0);
V_Svar_ts = np.mean(VImg_Svar_ts,axis=0);
V_Evar_ts = np.mean(VImg_Evar_ts,axis=0);
####### Whole Variances

V_w_Avar  = np.sum(V_Avar_ts);
V_w_Dvar  = np.sum(V_Dvar_ts);
V_w_Svar  = np.sum(V_Svar_ts);
V_w_Evar  = np.sum(V_Evar_ts);

###########Global
V_g_Avar_ts = bar_Y**2;
V_g_Dvar_ts = bar_Dvar**2/4;
V_g_Svar_ts = bar_Svar**2/4;
V_g_Evar_ts = bar_Y[np.array((0,T-1))]**2/2;

###########Global ts (Just for vis)
#V_g_Dts=bar_Dvar/2;
#V_g_Ats=bar_Y;
#V_g_Sts=bar_Svar/2;

V_g_Avar  = np.sum(V_g_Avar_ts);
V_g_Dvar  = np.sum(V_g_Dvar_ts);
V_g_Svar  = np.sum(V_g_Svar_ts);
V_g_Evar  = np.sum(V_g_Evar_ts);

#V.g_Avar  = np.sum(B.Ybar**2);
#V.g_Dvar  = np.sum(B.Dbar**2)/4;
#V.g_Svar  = np.sum(B.Sbar**2)/4;
#V.g_Evar  = np.sum(B.Ybar([1,T0])**2)/2;

#sys.exit("Stop HERE")
#exit()
#Here Here #Here Here #Here Here
###########Non-Global
V_ng_Avar_ts = np.mean((Y-np.tile(bar_Y,[I1,1]))**2,axis=0);
V_ng_Dvar_ts = np.mean((Dvar-np.tile(bar_Dvar,[I1,1]))**2,axis=0)/4;
V_ng_Svar_ts = np.mean((Svar-np.tile(bar_Svar,[I1,1]))**2,axis=0)/4;
V_ng_Evar_ts = np.mean([(Yhead-bar_Y1bar)**2,(Ytail-bar_Ytbar)**2],axis=0)/2;

V_ng_Avar = np.sum(V_ng_Avar_ts);
V_ng_Dvar = np.sum(V_ng_Dvar_ts);
V_ng_Svar = np.sum(V_ng_Svar_ts);
V_ng_Evar = np.sum(V_ng_Evar_ts);

vts = {'Avar_ts':V_Avar_ts,'Dvar_ts':V_Dvar_ts,
       'Svar_ts':V_Svar_ts,'Evar_ts':V_Evar_ts};

##SAVE Variance Time Series    
SaveMe2txt(vts,NewDir2DSE+"/VarTS")
##
    
vv = {'Avar':V_w_Avar,'Dvar':V_w_Dvar,'Svar':V_w_Svar,'Evar':V_w_Evar,
       'gAvar': V_g_Avar ,'gDvar': V_g_Dvar,'gSvar': V_g_Svar,'gEvar': V_g_Evar,
       'ngAvar':V_ng_Avar,'ngDvar':V_ng_Dvar,'ngSvar':V_ng_Svar,'ngEvar':V_ng_Evar};       

print "+DSE Variances:"
for keys, values in vv.items():
        print keys + "      : " + str(values)
        print "_________________________________"

         
#V.GrandMean_Untouched  = np.mean(mvY_Untouched);
#V.GrandMean_NormInt    = np.mean(mvY_NormInt);
#V.GrandMean_Demeaned   = np.mean(mvY_Demeaned);
#V.GranMean_WholeBrain  = np.mean(mvY_WholeImage);

#V.ng_Avar = sum(mean((Y-repmat(B.Ybar,[I1,1])).^2));
#V.ng_Dvar = sum(mean((D-repmat(B.Dbar,[I1,1])).^2))./4;
#V.ng_Svar = sum(mean((S-repmat(B.Sbar,[I1,1])).^2))./4;
#V.ng_Evar = sum(mean([(Yhead-B.Y1bar).^2,(Ytail-B.Ytbar).^2]))./2;

# Sanity Chek - The moment of truth!
# gvars             = V.g_Dvar+V.g_Svar+V.g_Evar;
# rgvars            = V.rg_Dvar+V.rg_Svar+V.rg_Evar;
# WholeWholeVar_Test=gvars+rgvars;
# %assert(WholeWholeVar==WholeWholeVar_Test,'VarDecomp failed')
# disp(['WholeVar= ' num2str(V.w_Avar) ' , sum of decomp var= ' num2str(WholeWholeVar_Test)])

########### SED ANOVA table
SS      = I1 * np.array((V_w_Avar,V_w_Dvar,V_w_Svar,V_w_Evar,V_g_Avar,V_g_Dvar,V_g_Svar,V_g_Evar))
MS      = SS/I1/T
RMS     = np.sqrt(MS)
Prntg   = RMS**2/RMS[0]**2*100
Expd    = np.concatenate((np.array((1,(T-1)/T/2,(T-1)/T/2,1/T),dtype='float64'),np.array((1,(T-1)/T/2,(T-1)/T/2,1/T),dtype='float64')/I1))
RelVar  = Prntg/100/Expd

SaveMe2txt({'Percents':Prntg,'Relatives':RelVar,'SumOfSquares':SS,'RootMeanOfSquares':MS},NewDir2DSE+"/VarStat")

DSEVarfig = plt.figure()       
plt.plot(vts['Svar_ts'])
plt.plot(vts['Dvar_ts'])
plt.legend(['Svar_ts', 'Dvar_ts'], loc='upper left')
plt.savefig(NewDir2DSE+"/VarStat/DSvars.png")
plt.close(DSEVarfig)
#plt.show() 

#Var_Tab = [V_w_Avar,V_w_Dvar,V_w_Svar,V_w_Evar;...
#    V_g_Avar,V_g_Dvar,V_g_Svar,V_g_Evar;...
#    V_ng_Avar,V_ng_Dvar,V_ng_Svar,V_ng_Evar];

###############################################################################

