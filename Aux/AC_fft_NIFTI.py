#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 16:53:44 2018

@author: sorooshafyouni
"""

def AC_fft_NIFTI(Y)
    Y = Y-np.transpose(np.tile(mY2,(T,1)));  #deamean to be safe!

    nfft    = 2**nextpow2(2*L-1); #zero-pad the hell out!
    yfft    = np.fft.fft(Y,n=nfft,axis=1); #be careful with the dimensions
    ACOV = np.real(np.fft.ifft(yfft*np.conj(yfft),axis=1));
    ACOV = ACOV[:,1:L];
    xAC  = ACOV/np.sum(np.abs(Y)**2,axis=1); #normalise the COVs   
    #bnd=(sqrt(2)*erfinv(0.95))./sqrt(L); #assumes normality for AC
    #CI=[-bnd bnd];
    return(xAC)