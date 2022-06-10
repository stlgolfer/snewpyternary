#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 23:16:27 2022
This contains data handlers for handling individual model data
@author: phyics
"""

import numpy as np
import os

supported_detectors = [
    'scint20kt',
    'ar40kt',
    'wc100kt30prct'
    ]

def build_detector_profiles():
    '''
    Builds a profile for a detector that includes a dictionary of data handlers
    and the corresponding data labels (used for labeling axes)

    Returns
    -------
    detector_defs : dictionary
        The dictionary of profiles

    '''
    detector_defs = {}
    module_scope = __import__(__name__)
    for detector in supported_detectors:
        detector_defs[detector] = {
            'handler':getattr(module_scope, f'h_{detector}'),
            'axes':getattr(module_scope,f'axes_{detector}')
            }
    return detector_defs

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

'''
Data handler axes labels definitions
'''

def axes_scint20kt():
    return ['ibd','nue+es','nc']

def axes_ar40kt():
    return ['nue','nuebar','nc']

def axes_wc100kt30prct():
    return ['ibd','nue+es','nc']