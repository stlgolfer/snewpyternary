#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 23:16:27 2022
This contains data handlers for handling individual model data
@author: phyics
"""

import numpy as np

def h_scint20kt(data):
    # must return a list of a, b, c
    ibd = np.sum(data['ibd'])
    nue_plus_es=np.sum(data['nue_C12'])+np.sum(data['nue_C13']+data['e'])
    nc = np.sum(data['nc'])
    return [ibd,nue_plus_es,nc]

def h_ar40kt(data):
    return [np.sum(data['nue_Ar40']),np.sum(data['nuebar_Ar40']),np.sum(data['nc'])]

def h_wc100kt30prct(data):
    ibd = np.sum(data['ibd'])
    nue_plus_es=np.sum(data['nue_O16'])+np.sum(data['e'])
    nc = np.sum(data['nc'])
    return [ibd,nue_plus_es,nc]