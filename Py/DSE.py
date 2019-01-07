#!/usr/bin/env python3
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
    
def CleanNIFTI(V,\
               scl = 1,\
               demean = True,\
               **kwargs):
    
    import numpy as np 
    import nibabel as nib #pip install nibabel
    import os
    
    if isinstance(V, str):
        print('Input is a path to a nifti file')
        [pth, fn] = os.path.split(V)
        
        filename =  os.path.expandvars(V)
        imobj = nib.load(filename, mmap=False)
        imdat = imobj.get_data().astype(float)
        
        T = imdat.shape[3]
        I = np.prod(imdat.shape[0:3])
        Y = np.reshape(imdat,(I,T))
                   
    else:
        print('Input is a matrix.')
        Y = V; del V
        T = np.shape(Y)[1]
        I = np.shape(Y)[0]
        imobj = ''


    print("The image has total of: " + str(I) + " voxels & " + str(T) + " time-points.")
    #find zeros or nans
    zrIdx = np.where(~Y.any(axis=1))[0];
    nanIdx = np.where(np.isnan(Y))[0]; 
    rmvIdx = np.concatenate((zrIdx,nanIdx));
    
    Y = np.delete(Y,rmvIdx,axis=0)
    I1 = Y.shape[0]
    
    print("After cleaning; " + str(I1) + " voxels & " + str(T) + " time-points.")     
        
    #Intensity Normalisation############################## 
    print("Intensity Normalisation to " + str(scl))
    
    md = np.median(np.mean(Y,axis = 1)) #median of mean image
    Y = Y/md * scl;
    ######################################################
    
    #Demean each time series##############################
    if demean: 
        print("Demean along T")
        mY2 = np.mean(Y,axis=1)
        Y = Y-np.transpose(np.tile(mY2,(T,1)));
    ######################################################            
    
        
    return (Y, T, I, I1, rmvIdx, imobj)

######################################################
######################################################
def SaveMe2txt(Dict2Save,Where2Save):
    import numpy as np 
    import os
    
    if not os.path.exists(Where2Save):
        os.makedirs(Where2Save)
    for varname in Dict2Save.keys():       
        np.savetxt(Where2Save + '/DSE_' + str(varname) + '.txt',Dict2Save[varname],fmt='%10.5f')
######################################################
######################################################
def SaveMe2Nifti(Dict2Save,Where2Save,imobj,rmvIdx):
    import numpy as np 
    import nibabel as nib #pip install nibabel
    import os
    
    if not os.path.exists(Where2Save):
        os.makedirs(Where2Save)

        
    imghdr = imobj.header.copy()
    [X,Y,Z,T] = imobj.shape[0:4]
    I = np.prod((X,Y,Z))
    idx_tmp  = range(I+1)[1:I+1]
    idx_tmp0 = np.delete(idx_tmp,rmvIdx,axis=0)        
    
    for varname in Dict2Save.keys():    
        print('writing *' + varname + '* images...')
        
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
            
        else: raise(ValueError ,"Something is wrong with image dimensions.")
        
        imgobj_vars = nib.nifti1.Nifti1Image(I3Y, None, header=imghdr)
        nib.save(imgobj_vars, Where2Save+'/DSE_'+ str(varname) +'.nii.gz')
        
        del(Y2I,I3Y,G) 
##########################################################
##########################################################

def DSE_Calc(Y,\
        NewDir2DSE='',\
        scl = 1,\
        demean = True,\
        plotme = True,\
        **kwargs):
    import numpy as np 
    import matplotlib.pyplot as plt
    import os
    
    if NewDir2DSE=='':
        NewDir2DSE = (os.getcwd() + '/DSE/')
        
    #READ AND CLEAN###########################################
    
    NFT  = CleanNIFTI(Y,scl=scl,demean=demean)
    Y = NFT[0];
    T = NFT[1];
    #I = NFT[2];
    I1 = NFT[3];
    rmvIdx = NFT[4];
    imobj  = NFT[5];
    ##########################################################
    
    
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
    print("Variance 4D Images")
    VImg_Avar_ts = Y**2;
    VImg_Dvar_ts = Dvar**2/4;
    VImg_Svar_ts = Svar**2/4;
    VImg_Evar_ts = np.transpose(np.array((Yhead,Ytail))**2/2);
    
    #Averaged across time#######################################
    print("Variance 3D Images")
    VImg_Avar = np.mean(Y**2,axis=1);
    VImg_Dvar = np.mean(Dvar**2,axis=1)/2;
    VImg_Svar = np.mean(Svar**2,axis=1)/2;
    VImg_Evar = np.mean(np.array((Yhead**2,Ytail**2)),axis=0); # <<<< should be checked
    
    if imobj != '':
        print('Writing the DSE images...')
        vimg = {'Avar_Img3':VImg_Avar,'Dvar_Img3':VImg_Dvar,
            'Svar_Img3':VImg_Svar,'Evar_Img3':VImg_Evar,
            'Avar_Img4':VImg_Avar_ts,'Dvar_Img4':VImg_Dvar_ts,
            'Svar_Img4':VImg_Svar_ts,'Evar_Img4':VImg_Evar_ts};
    
        SaveMe2Nifti(vimg,NewDir2DSE+"/VarImg",imobj,rmvIdx)
    else: vimg=''
             
    #DSE Time series -- averaged across I#######################
    print("DSE time series")
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
    
    print("+DSE Variances:")
    for keys, values in vv.items():
            print (keys + "      : " + str(values))
            print ("_________________________________")
    
             
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
    SS      = I1 * np.array((V_w_Avar,V_w_Dvar,V_w_Svar,V_w_Evar,\
                             V_g_Avar,V_g_Dvar,V_g_Svar,V_g_Evar))
    MS      = SS/I1/T
    RMS     = np.sqrt(MS)
    Prntg   = RMS**2/RMS[0]**2*100
    Expd    = np.concatenate(\
                             (np.array((1,(T-1)/T/2,(T-1)/T/2,1/T),dtype='float64'),\
                              np.array((1,(T-1)/T/2,(T-1)/T/2,1/T),dtype='float64')/I1))
    RelVar  = Prntg/100/Expd
    
    SaveMe2txt({'Percents':Prntg,\
                'Relatives':RelVar,\
                'SumOfSquares':SS,\
                'RootMeanOfSquares':MS},\
    NewDir2DSE+"/VarStat")
    
    if plotme:
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
    
    ######################################################
    return vv,vts,Prntg,MS,RMS,RelVar,Expd,vimg
##########################################################
##########################################################   
def DVARS_Calc(Y,\
               dd=1, \
               WhichExpVal = 'median', WhichVar = 'hIQRd',\
               scl = 1,\
               demean = True,\
               **kwargs): 
    import scipy.stats as st
    import numpy as np 
    #defaults:
    #WhichExpVal = 'median';    
    #WhichVar    = 'hIQRd';
    #dd = 1;
    
    ##FUNCS####################################################
    def IQRsd(x):   return ((np.percentile(x,75,axis=0)  - np.percentile(x,25,axis=0))/1.349);
    def H_IQRsd(x): return ((np.percentile(x,50,axis=0)  - np.percentile(x,25,axis=0))/1.349*2);
    #--
    #if tsflag
    #    def Zstat(x,m,s): abs((x-m)/s);
    #else
    
    def ZStat(x,m,s):   return((x-m)/s)
    def Zp(x,m,s):      return(1-st.norm.cdf(ZStat(x,m,s)))
    #--
    def X2Stat(x,m,s):  return(2*m/s**2*x)
    def X2df(m,s):      return(2*m**2/s**2)
    #X2p    =   @(x,m,s) 1-chi2cdf(X2stat(x,m,s),X2df(m,s));
    def X2p0(x,m,s):    return((X2Stat(x,m,s)-X2df(m,s))/np.sqrt(2*X2df(m,s)));
    def X2p(x,m,s):     return(st.chi2.cdf(X2Stat(x,m,s),X2df(m,s)));
    
    
    #READ AND CLEAN###########################################
    
    NFT  = CleanNIFTI(Y,scl=scl,demean=demean)
    Y = NFT[0];
    T = NFT[1];
    #I = NFT[2];
    I1 = NFT[3];
    
    ##########################################################
    
    #Dvar
    Dvar = Y[:,0:T-1]-Y[:,1:T] #OMG, python!
    
    DVARS = np.sqrt(np.sum(Dvar**2,axis=0)/I1)
    
    DVARS2  = np.mean(Dvar**2,axis=0);
    Rob_S_D = IQRsd(Dvar);
    
    Z       = DVARS2**dd;
    M_Z     = np.median(Z,axis=0);
    
    
    print('Estimating the moments...')
    DVMean = {'sig2bar':np.mean(Rob_S_D**2,axis=0),\
              'sig2median':np.median(Rob_S_D**2,axis=0),\
              'median':np.median(DVARS2,axis=0),\
              'sigbar2':np.mean(Rob_S_D,axis=0)**2,\
              'xbar':np.mean(DVARS2,axis=0)};   
              
    DVVar  = {'S2':np.var(DVARS2,axis=0),\
              'IQRd':(1/dd*M_Z**(1/dd-1)*IQRsd(Z))**2,\
              'hIQRd':(1/dd*M_Z**(1/dd-1)*H_IQRsd(Z))**2};

    print('Estimating Delta % Vars')
    DeltapDvar = (DVARS2-np.median(DVARS2,axis=0))/(4*np.mean(Y**2))*100; #This is one other Standardised DVARS!
    
    
    print('Doing the inference...')
    M_DV2         = DVMean[WhichExpVal];             
    S_DV2         = np.sqrt(DVVar[WhichVar]);       
    Pval          = X2p(DVARS2,M_DV2,S_DV2);
    Zval          = ZStat(DVARS2,M_DV2,S_DV2);
    c             = X2Stat(DVARS2,M_DV2,S_DV2);
    nu            = X2df(M_DV2,S_DV2); #Spatial EDF 
    NDVARS_X2     = -st.norm.ppf(Pval) #-norminv(Pval);    
    NDVARS_X20    = X2p0(DVARS2,M_DV2,S_DV2);
    
    #only substitute the infs, the rest is better done in matlab.
    WhereIsInf = np.where(np.isinf(NDVARS_X2))
    NDVARS_X2[WhereIsInf] = NDVARS_X20[WhereIsInf]  #This is one version of the standardised DVARS!
    
    

##########################################################
##########################################################

def main():
        ################################Read and Clean
        V = '/Users/sorooshafyouni/Desktop/HCP\
        /135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/\
        rfMRI_REST1_LR_hp2000_clean.nii.gz'

        ###################################################################
        NewDir2DSE = pth + '/DSE/'

        #Do the DSE Analysis and Save the file and everything        
        DSE(Y,I,T,I1,imobj,pth,fn,rmvIdx,Dvar,NewDir2DSE)
        
        
if __name__ == "__main__":
    main()   
