#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:31:58 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

from DSE import DSE_Calc, DVARS_Calc, CleanNIFTI
import pickle


RS_fMRI_Dir    = '/Users/sorooshafyouni/Desktop/HCP/135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz'
DSEResults_Dir = '/Users/sorooshafyouni/Home/GitClone/DVARS/Python/'


#Write me to a pickle, quickly!
def Save2Pickle(Path2File,Vars):
    with open(Path2File,"wb") as f:
        pickle.dump(Vars, f)		

##++++++++++++++++++++++++++++++++++++
##+++++++++++++++++++CLEAN THE NIFTI+++
#
print('========Read and clean...')

print('========Running DSE...')
##+++++++++++++++++++DSE+++++
DSEOut  = DSE_Calc( RS_fMRI_Dir,\
        scl         = 0,\
        norm        = 0,\
        DSEImage    = False,\
        verbose     = False)

#Vout = CleanNIFTI(RS_fMRI_Dir,scl = 1/100)
##
#GS = Vout[6] # this is the global signal 
#YY = Vout[0] # This is the cleaned version of the image. 
#
#del Vout
#
## Just to test everything is sound with the DSE and DAVRS code. 
##YY = np.random.randn(1500,500)
#
#
#
#print('========Running DSE...')
###+++++++++++++++++++DSE+++++
#DSEOut  = DSE_Calc( YY,\
#        scl         = 0,\
#        norm        = 0,\
#        DSEImage    = False,\
#        verbose     = False)
#
#
#
#
#print('========Running DVARS...')
##+++++++++++++++++++DVARS++++
#DVARSOut = DVARS_Calc(YY,\
#               dd            = 1, \
#               WhichExpVal   = 'median',\
#               WhichVar      = 'hIQRd',\
#               norm          = 0,\
#               scl           = 0,\
#               DeltapDvarThr = 5)
#
#
##++++Write me on a pickle!++++++++++
#Save2Pickle(DSEResults_Dir + 'GS_Demo.pkl',GS)
#Save2Pickle(DSEResults_Dir + 'DSE_Demo.pkl',DSEOut)
#Save2Pickle(DSEResults_Dir + 'DVARS_Demo.pkl',DVARSOut)
#
#print('DONE!')
