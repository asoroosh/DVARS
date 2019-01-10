#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 14:50:41 2019

@author: sorooshafyouni
University of Oxford, 2019
"""

def myfunc(V, foo="default value", **kwargs):
    #foo = kwargs.pop('foo')
    
    print(foo)
    print(V)