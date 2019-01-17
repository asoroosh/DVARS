#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:31:58 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

from DSE import DSE_Calc, DVARS_Calc, CleanNIFTI


V = '/Users/sorooshafyouni/Desktop/HCP/135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz'

#+++++++++++++++++++CLEAN THE NIFTI+++

Vout = CleanNIFTI(V,scl = 1/100)

#+++++++++++++++++++DSE+++

DSEOut  = DSE_Calc(Vout[0],\
        scl        = 0,\
        demean     = False,\
        DSEImage   = False,\
        verbose    = True)

#+++++++++++++++++++DVARS+++
DVARSOut = DVARS_Calc(Vout[0],\
               dd=1, \
               WhichExpVal = 'median', WhichVar = 'hIQRd',\
               scl = 0,\
               demean = False,\
               DeltapDvarThr = 5)