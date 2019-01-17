#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:31:58 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

from DSE import DSE_Calc, DVARS_Calc, CleanNIFTI
import numpy as np

V = '/Users/sorooshafyouni/Desktop/HCP/135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz'
#
##++++++++++++++++++++++++++++++++++++
##+++++++++++++++++++CLEAN THE NIFTI+++
#
Vout = CleanNIFTI(V,scl = 1/100)
#
YY = Vout[0]

#YY = np.random.randn(1500,500)

##+++++++++++++++++++++++++++
##+++++++++++++++++++DSE+++++
DSEOut  = DSE_Calc(YY,\
        scl        = 0,\
        norm       = 0,\
        DSEImage   = False,\
        verbose    = True)

#++++++++++++++++++++++++++++
#+++++++++++++++++++DVARS++++
DVARSOut = DVARS_Calc(YY,\
               dd = 1, \
               WhichExpVal  = 'median',\
               WhichVar     = 'hIQRd',\
               norm = 0,\
               scl  = 0,\
               DeltapDvarThr = 5)