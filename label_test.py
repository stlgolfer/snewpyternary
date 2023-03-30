#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 14:55:14 2022

@author: physics
"""

import snewpyternary as t
from snewpy.models import Nakazato_2013
from astropy import units as u

modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"

t.create_default_detector_plot([(40,50,10)], ['top','right','left'], 'Plot test')
t.create_default_detector_plot([(1,1,98)], ['top','right','left'], 'Plot test 2')
t.create_default_detector_plot([(1,98,1)], ['top','right','left'], 'Plot test 2')
t.create_default_detector_plot([(98,1,1)], ['top','right','left'], 'Plot test 2')

# p_data, r_data = t.create_flux_scatter(modelFilePath, model_type, model, 0.1*u.s, log_bins=True, use_cache=True)
#
# t.create_regular_plot(r_data, ['nux','anue','nue'], 'Regular plot test', 'Numbers')
# t.create_default_flux_plot(p_data, 'Regular plot test')