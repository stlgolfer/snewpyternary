#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 23:43:12 2022

@author: physics
"""

from snewpy.models import Nakazato_2013, Bollig_2016

class SNEWPYModel:
    """
    SNEWPY models have several available configs, so it's convenient to just
    bundle this information togeter
    """
    
    def __init__(self,file_path,model_type,model):
        self.file_path = file_path
        self.model_type = model_type
        self.model = model

# now we must process each one

snewpy_models = {}

'''
Nakazato_2013
'''
snewpy_models['Nakazato_2013'] = SNEWPYModel(
    './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits',
    'Nakazato_2013',
    Nakazato_2013('./SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits')
    )

'''
Bollig_2016
'''
fname = './SNEWPY_models/Bollig_2016/s11.2c'
snewpy_models['Bollig_2016'] = SNEWPYModel(
    fname,
    'Bollig_2016',
    Bollig_2016(fname)
    )

'''
Walk_2019
'''
fname = './SNEWPY_models/Walk_2019/s40.0c_3DBH_dir1'
snewpy_models['Walk_2019'] = SNEWPYModel(
    fname,
    'Walk_2019',
    Bollig_2016(fname)
    )