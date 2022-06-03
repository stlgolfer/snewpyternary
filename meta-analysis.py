#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 01:45:12 2022
Here we aim to use the parametrizations of the library to create example
plots from the SNEWPY data
@author: phyics
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.models import Nakazato_2013
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
import snewpyternary as t
import os
import ternary
import math
from snewpy.flavor_transformation import *

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

#simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
deltat = (1*u.s) # time bin size, details found in SNEWPY article
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = True
weighting = "unweighted"
model

detector_list = ['scint20kt','ar40kt','wc100kt30prct']
# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
transform_list = list(flavor_transformation_dict.keys())

# for d in detector_list:
#     # for each detector config, simulate all types of transformations
#     for transform in transform_list:
#         plot_title = "Nakazato_2013 Events " + weighting + " "  + d + " " + transform
#         plot_data = t.create_detector_event_scatter(modelFilePath,model_type, d, model,transformation=transform)
#         # then we can generate a ternary plot for this
#         t.create_default_detector_plot(plot_data,plot_title)

# create flux plot as well
for transform in transform_list:
    flux_scatter_data = t.create_flux_scatter(modelFilePath, model_type, model,transform=transform)
    t.create_default_flux_plot(flux_scatter_data, "{model} Flux {transform}".format(model=model_type,transform=transform))