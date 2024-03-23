#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DSE decomposition and DVARS inference. 

REFERENCE:
    Afyouni, Soroosh, and Thomas E. Nichols. "Insight and Inference for DVARS."
    NeuroImage 172 (2018): 291-312.

REQUIREMENT:
    os, numpy, nibabel, scipy
    
FUNCTIONS:
    CleanNIFTI, DSE_Calc, DVARS_Calc
    
Created on Mon Mar 12 11:56:21 2018
@author: sorooshafyouni
University of Oxford, 2018
srafyouni@gmail.com
"""
def CleanNIFTI(V,\
               scl  = 0,\
               norm = 0,\
               demean = True,\
               copy   = True):
    """
    CleanNIFTI(V,scl = 1,demean = True)
    Cleans and reshape the input images. 
    
    INPUT:
        V: Can be (1) a string indicating the path to the 
        nifti file (2) a numerical matrix of size IxT. 
        Where I is number of voxels (I=Nx x Ny x Nz) and T is 
        number of data-points.
           
        scl [optional]: Scaling factor in intensity normalisation. To down 
        scale us fractions. [default: 0, meaning there will be 
        no scaling]
        
        demean [optional]: if True it will demean the time series. [default: True]
    
    SA, Ox, 2019
    srafyouni@gmail.com
    """
    import numpy as np
    import nibabel as nib #pip install nibabel
    import os
    
    if copy:
        #Also V is going to be copied only if it is not str
        #see the else argument in the next if statement
        scl  =  np.copy(scl)
        norm =  np.copy(norm)

    if isinstance(V, str):
        print('CleanNIFTI::: Input is a path to a nifti file')
        [pth, fn] = os.path.split(V)

        filename =  os.path.expandvars(V)
        imobj = nib.load(filename, mmap=False)
        imdat = imobj.get_data().astype(float)

        T = imdat.shape[3]
        I = np.prod(imdat.shape[0:3])
        Y = np.reshape(imdat,(I,T))

    else:
        Y = V.copy()
        del V #to save some memory, ffs
        shapes = np.shape(Y)
        if len(shapes)==2:
            print('CleanNIFTI::: Input is a 2D matrix.')
            T = np.shape(Y)[1]
            I = np.shape(Y)[0]
            imobj = ''
        elif len(shapes)==4:
            print('CleanNIFTI::: Input is a 4D matrix.')
            I = np.prod(np.shape(Y)[:3])
            T = np.shape(Y)[-1]
            Y = np.reshape(Y,(I,T))
        else: raise ValueError('CleanNIFTI::: the shape of the input should be either 2D or 4D')

    print("CleanNIFTI::: The image has total of: " + str(I) + " voxels & " + str(T) + " time-points.")
    #find zeros or nans
    zrIdx  = np.where(~Y.any(axis=1))[0];
    nanIdx = np.where(np.isnan(Y))[0];
    rmvIdx = np.concatenate((zrIdx,nanIdx));

    Y = np.delete(Y,rmvIdx,axis=0)
    I1 = Y.shape[0]

    print("CleanNIFTI::: After cleaning; " + str(I1) + " voxels & " + str(T) + " time-points.")

    #Intensity Normalisation##############################
    
    MeanImage = np.mean(Y,axis = 1)
    if scl!=0 and norm==0:
        print("CleanNIFTI::: Scaled by " + str(scl))
        #md = np.median(MeanImage) #median of mean image
        md = 1
        Y  = Y/md * scl;        
    elif norm!=0 and scl==0:
        print("CleanNIFTI::: Intensity Normalisation done by: " + str(scl))
        md = np.median(MeanImage) #median of mean image
        Y  = Y/md * norm;
    elif norm==0 and scl==0: 
        print('CleanNIFTI::: No normalisation/scaling has been set!')
    else:
        raise ValueError('CleanNIFTI::: norm (' + str(norm) + ') and scl (' + str(scl) + ') parameters are wrong!')
        
    #Global Signal#########################################
    GS = np.mean(Y,axis = 0)        
        
    #Demean each time series##############################
    if demean:
        print("CleanNIFTI::: Demean along T")
        Y = Y-np.transpose(np.tile(MeanImage,(T,1)));

#    return {'Y':Y,\
#            'T':T,\
#            'I':I,\
#            'I1':I1,\
#            'removables':rmvIdx,\
#            'ImageObj':imobj}

    return (Y,T,I,I1,rmvIdx,imobj,GS,MeanImage)

######################################################
######################################################
def SaveMe2txt(Dict2Save,Where2Save):
    """ Save the input dictionary to a text file. SA, OX, 2019"""
    import numpy as np
    import os

    if not os.path.exists(Where2Save):
        os.makedirs(Where2Save)
        
    for varname in Dict2Save.keys():
        np.savetxt(Where2Save + '/DSE_' + str(varname) + '.txt',Dict2Save[varname],fmt='%10.5f')
######################################################
######################################################
def SaveMe2Nifti(Dict2Save,Where2Save,imobj,rmvIdx):
    """ Save the input dictionary (VarImages) to a NIFTI files.
    
    SA, OX, 2019"""
    
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
        print('SaveMe2Nifti::: writing *' + varname + '* images...')

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

        else: raise(ValueError ,"SaveMe2Nifti::: Something is wrong with image dimensions.")

        imgobj_vars = nib.nifti1.Nifti1Image(I3Y, None, header=imghdr)
        nib.save(imgobj_vars, Where2Save+'/DSE_'+ str(varname) +'.nii.gz')

        del(Y2I,I3Y,G)
##########################################################
##########################################################


def DSE_Calc(Y,\
        NewDir2DSE = '',\
        scl        = 0,\
        norm       = 0,\
        demean     = True,\
        DSEImage   = False,\
        imobj      = '',\
        verbose    = False):
    """
    DSE_Calc(Y, NewDir2DSE='', scl = 0, demean = True, \
             DSEImage = False, verbose = False)    
    Calculate the DSE descomposition components; fast slow and edge variabilities, 
    in form of scalar, time series, 3D and 4D images. 

INPUTS:    
    
    V0:            Can be (1) a string indicating the path to the 
                   nifti/cifti file (2) a numerical matrix of size IxT. 
                   Where I is number of voxels (I=Nx x Ny x Nz) and T is 
                   number of data-points.

   'scl':          Output directory. Should only be used when the
                   input is a nifti image and user needs to save the
                   S, D and E (3D and 4D) images. [this is a param related to 
                   function CleanNIFTI]
                   
    demean [optional]: if True it will demean the time series. [default: True]               
                   
    'NewDir2DSE'   Path to directory which contain the results; including tables 
                   and variance images
                   e.g.: DSEOut=dse.DSE_Calc(V0,'NewDir2DSE'='~/Where/to/save/')
                   
   'DSEImage'      If set True then the function will reconstruct 3D and 4D 
                   images and save them. NB! It requires lots of storage 
                   [defualt: Flase]
   
   'verbose'     : Set to 1 if you need the function to save ANOVA tables, plots
                   and summary files.[default:False]
    
    Soroosh Afyouni, Ox, 2019
    srafyouni@gmail.com
    """
    import numpy as np
    import os

    if NewDir2DSE=='':
        NewDir2DSE = (os.getcwd() + '/DSE/')

    #READ AND CLEAN###########################################
    
    NFT  = CleanNIFTI(Y,scl=scl,demean=demean,norm=norm)
    Y = NFT[0];
    T = NFT[1];
    #I = NFT[2];
    I1     = NFT[3];
    rmvIdx = NFT[4];
    imobj  = NFT[5];
    ##########################################################

    bar_Y  = np.sum(Y,axis=0)/I1; #global signal is here!

    #Dvar
    Dvar = Y[:,0:T-1]-Y[:,1:T] #OMG, python!
    bar_Dvar  = np.sum(Dvar,axis=0)/I1;

    #print(Dvar[:5,:5])
    #print(I1)

    #Svar
    Svar = Y[:,0:T-1]+Y[:,1:T]
    bar_Svar = np.sum(Svar,axis=0)/I1;

    #print(Svar[:5,:5])
    #print(I1)

    #Evar
    Ytail = Y[:,-1];
    Yhead = Y[:,0];
    bar_Ytbar = np.sum(Ytail)/I1;
    bar_Y1bar = np.sum(Yhead)/I1;

    #SED Var Images######################################################
    print("DSE_Calc::: Variance 4D Images")
    VImg_Avar_ts = Y**2;
    VImg_Dvar_ts = Dvar**2/4;
    VImg_Svar_ts = Svar**2/4;
    VImg_Evar_ts = np.transpose(np.array((Yhead,Ytail))**2/2);

    #Averaged across time#######################################
    print("DSE_Calc::: Variance 3D Images")
    VImg_Avar = np.mean(Y**2,axis=1);
    VImg_Dvar = np.mean(Dvar**2,axis=1)/2;
    VImg_Svar = np.mean(Svar**2,axis=1)/2;
    VImg_Evar = np.mean(np.array((Yhead**2,Ytail**2)),axis=0); # <<<< should be checked

    if imobj != '' and DSEImage:
        print('DSE_Calc::: Writing the DSE images...')
        vimg = {'Avar_Img3': VImg_Avar,\
                'Dvar_Img3': VImg_Dvar/VImg_Avar*100,\
                'Svar_Img3': VImg_Svar/VImg_Avar*100,\
                'Evar_Img3': VImg_Evar/VImg_Avar*100,\
                'Avar_Img4': VImg_Avar_ts,\
                'Dvar_Img4': VImg_Dvar_ts/np.transpose(np.tile(VImg_Avar,(T-1,1))*100),\
                'Svar_Img4': VImg_Svar_ts/np.transpose(np.tile(VImg_Avar,(T-1,1))*100),\
                'Evar_Img4': VImg_Evar_ts/np.transpose(np.tile(VImg_Avar,(2,1))*100)};

        SaveMe2Nifti(vimg,NewDir2DSE+"/VarImg",imobj,rmvIdx)
    else:
        vimg=''
        print('DSE_Calc::: The DSE image option was left False or there is no image object, if you need them trigger option DSEImage=.')

    #DSE Time series -- averaged across I#######################
    print("DSE_Calc::: DSE time series")
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

    vts = {'Avar_ts':V_Avar_ts,\
           'Dvar_ts':V_Dvar_ts,\
           'Svar_ts':V_Svar_ts,\
           'Evar_ts':V_Evar_ts};

    normvts = {'pDvar_ts':V_Dvar_ts/np.mean(V_Avar_ts),\
               'pSvar_ts':V_Svar_ts/np.mean(V_Avar_ts),\
               'pEvar_ts':V_Evar_ts/np.mean(V_Avar_ts)};


    vv = {'Avar'  : V_w_Avar,\
          'Dvar'  : V_w_Dvar,\
          'Svar'  : V_w_Svar,\
          'Evar'  : V_w_Evar,\
          'gAvar' : V_g_Avar ,\
          'gDvar' : V_g_Dvar,\
          'gSvar' : V_g_Svar,\
          'gEvar' : V_g_Evar,\
          'ngAvar': V_ng_Avar,\
          'ngDvar': V_ng_Dvar,\
          'ngSvar': V_ng_Svar,\
          'ngEvar': V_ng_Evar};

    print("DSE_Calc::: DSE Variances:")
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



    if verbose:
        import matplotlib.pyplot as plt
        
        VarStatDir = NewDir2DSE+'/VarTS'
        if not os.path.exists(VarStatDir):
            os.makedirs(VarStatDir)
        
        DSEVarfig = plt.figure()
        plt.plot(vts['Svar_ts'])
        plt.plot(vts['Dvar_ts'])
        plt.legend(['Svar_ts', 'Dvar_ts'], loc='upper left')
        plt.savefig(VarStatDir+"/DSvars.png")
        plt.close(DSEVarfig)
        #plt.show()
        
        SaveMe2txt({'Percents':Prntg,\
                    'Relatives':RelVar,\
                    'SumOfSquares':SS,\
                    'RootMeanOfSquares':MS},\
        VarStatDir)        
        
        ##SAVE Variance Time Series
        SaveMe2txt(vts,VarStatDir)
        ##

        #Var_Tab = [V_w_Avar,V_w_Dvar,V_w_Svar,V_w_Evar;...
        #    V_g_Avar,V_g_Dvar,V_g_Svar,V_g_Evar;...
        #    V_ng_Avar,V_ng_Dvar,V_ng_Svar,V_ng_Evar];

    ######################################################
    DSEOut = {'VarScalar':vv\
            ,'VarTimeSeries':vts\
            ,'NormVarTimeSeries':normvts\
            ,'percent':Prntg\
            ,'MeanSquared':MS\
            ,'RMS':RMS\
            ,'Relative':RelVar\
            ,'Expected':Expd\
            ,'VarImages':vimg}
    
    return DSEOut

####################################################################################################################
####################################################################################################################

def DVARS_Calc(Y,\
               dd=1, \
               WhichExpVal = 'median', 
               WhichVar = 'hIQRd',\
               scl  = 0,\
               norm = 0,\
               demean = True,\
               DeltapDvarThr = 5,\
               copy = True):
    
    """
    DVARS_Calc(Y, dd=1, WhichExpVal = 'median', WhichVar = 'hIQRd',\
               scl = 0,demean = True,DeltapDvarThr = 5)    
    
    Find p-values, and performs statistical and practical significance testing 
    to identify corrupted volumes in a scan. Also calculates standardised DVARS

INPUTS:    
    
    V0:             Can be (1) a string indicating the path to the 
                    nifti/cifti file (2) a numerical matrix of size IxT. 
                    Where I is number of voxels (I=Nx x Ny x Nz) and T is 
                    number of data-points.
                   
   'dd':            Power of transformation [default:1]

    'WhichExpVal':  Method for robust estimate of expected value. The value
                    should be a digit corresponding to the order of 
                    following methods [default:'median']: 
                    'sig2bar','sig2median','median','sigbar2','xbar'.
                    For example: WhichExpVal='median' means the method 
                    to estimate robust expected value is empirical median. 
                    e.g. DVARS=dse.DVARS_Calc(V0,WhichExpVal='median')

    'WhichVar':    Method for robust estimate of variance. It can be either 
                   'IQR' for full IQR or 'hIQR' for half-IQR. 
                   [default:'hIQR']
                   e.g. DVARS=dse.DVARS_Calc(V0,WhichVar='hIQR')

    'DeltapDvarThr': Threshold (in percentage) for DeltapDvar variable which 
                     result in practical significance. [default: 5]      

   'scl':           Output directory. Should only be used when the
                    input is a nifti image and user needs to save the
                    S, D and E (3D and 4D) images. [this is a param related to 
                    function CleanNIFTI]
                   
   'demean':        if True it will demean the time series. [default: True]
   
   'verbose'     : Set to 1 if you need the function to save ANOVA tables, plots
                   and summary files.[default:False]
    
    Soroosh Afyouni, Ox, 2019
    srafyouni@gmail.com
    """
    
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

    def X2p0(x,m,s):    return((X2Stat(x,m,s)-X2df(m,s))/np.sqrt(2*X2df(m,s)));
    def X2p(x,m,s):     return(1-st.chi2.cdf(X2Stat(x,m,s),X2df(m,s)));


    #READ AND CLEAN###########################################

    NFT  = CleanNIFTI(Y,scl=scl,norm=norm,demean=demean)
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


    print('DVARS_Calc::: Estimating the moments...')
    DVMean = {'sig2bar':np.mean(Rob_S_D**2,axis=0),\
              'sig2median':np.median(Rob_S_D**2,axis=0),\
              'median':np.median(DVARS2,axis=0),\
              'sigbar2':np.mean(Rob_S_D,axis=0)**2,\
              'xbar':np.mean(DVARS2,axis=0)};

    DVVar  = {'S2':np.var(DVARS2,axis=0),\
              'IQRd':(1/dd*M_Z**(1/dd-1)*IQRsd(Z))**2,\
              'hIQRd':(1/dd*M_Z**(1/dd-1)*H_IQRsd(Z))**2};

    print('DVARS_Calc::: Estimating Delta % Vars')
    DeltapDvar = (DVARS2-np.median(DVARS2,axis=0))/(4*np.mean(Y**2))*100; #This is one other Standardised DVARS!


    print('DVARS_Calc::: Doing the inference...')
    M_DV2         = DVMean[WhichExpVal];
    S_DV2         = np.sqrt(DVVar[WhichVar]);
    
    Pval          = X2p(DVARS2,M_DV2,S_DV2);
    Zval          = ZStat(DVARS2,M_DV2,S_DV2);
    
    c             = X2Stat(DVARS2,M_DV2,S_DV2);
    nu            = X2df(M_DV2,S_DV2); #Spatial EDF
    
    NDVARS_X2     = -st.norm.ppf(Pval) # This should be checked!!
    NDVARS_X20    = X2p0(DVARS2,M_DV2,S_DV2);

    #only substitute the infs, the rest is better done in matlab.
    WhereIsInf = np.where(np.isinf(NDVARS_X2))
    NDVARS_X2[WhereIsInf] = NDVARS_X20[WhereIsInf]  #This is one version of the standardised DVARS!

    Alp_level = (0.05/(T-1));
    Hval_stat = np.where(Pval<Alp_level)
    print('DVARS_Calc::: ' + str(np.size(Hval_stat)) + ' of volumes are statistically significant.')

    Hval_practical = np.where(DeltapDvar>DeltapDvarThr)
    print('DVARS_Calc::: ' + str(np.size(Hval_practical)) + ' of volumes are practically significant.')

    Hval = np.intersect1d(Hval_stat,Hval_practical);
    print('DVARS_Calc::: ' + str(np.size(Hval)) + ' of volumes are significant.')

    #DVARSfig = plt.figure()
    #plt.plot(DeltapDvar)
    #plt.close(DVARSfig)
    #plt.show()

    DVARSOut = {'DVARS': {'DVARS': DVARS, 'DeltapDvar': DeltapDvar,'NDVARS_X2': NDVARS_X2},\
                'Inference': {'Pval': Pval,'H': Hval, 'HStat': Hval_stat, 'HPrac': Hval_practical},\
                'Stats': {'Mean': M_DV2,'SD':S_DV2,'DF':nu}}
    return DVARSOut

##########################################################
##########################################################
#def main():
#        ################################Read and Clean
#        V = '/Users/sorooshafyouni/Desktop/HCP\
#        /135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/\
#        rfMRI_REST1_LR_hp2000_clean.nii.gz'
#
#        ###################################################################
#        NewDir2DSE = pth + '/DSE/'
#
#        #Do the DSE Analysis and Save the file and everything
#        DSE(Y,I,T,I1,imobj,pth,fn,rmvIdx,Dvar,NewDir2DSE)
#
#
#if __name__ == "__main__":
#    main()
