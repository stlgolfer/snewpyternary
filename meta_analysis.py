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
from snewpy.models import Nakazato_2013, OConnor_2015
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
import snewpyternary as t
import os
import ternary
import math
from snewpy.flavor_transformation import *
import data_handlers as handlers
from multiprocessing import Process
import caching as cache

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

#simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
step = 0.01 #0.01 for best results
deltat=step*u.s
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'
multithreading = True # only use this when you don't need output
print(f'Timeframe of model is from {model.time[0]} to {model.time[len(model.time)-1]}')

# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
transform_list = list(flavor_transformation_dict.keys())

profiles = handlers.build_detector_profiles()

# let's do some parallel processsing to speed things up
def process_detector(detector, transform):
    plot_data, raw_data, l_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                                detector,
                                                model,
                                                deltat=deltat,
                                                transformation=transform,
                                                data_calc=profiles[detector]['handler'],
                                                use_cache=True
                                                )
    figure, tax = t.create_default_detector_plot(plot_data,
                                                  profiles[detector]['axes'](),
                                                  f'{model_type} {detector} {transform} Ternary',
                                                  save=True)
    # create left-point time bins
    time_bins = []
    for point in range(len(raw_data)):
        coordinate = (point*step)+0.5
        time_bins.append(coordinate)
        #if coordinate < 1: # this is for limiting the time domain
            
    t.create_regular_plot(
        plot_data=raw_data[:len(time_bins)],
        axes_titles=profiles[detector]['axes'](),
        plot_title=f'{model_type} {detector} {transform}',
        ylab="Event Counts",
        xlab="Right Time Bin Coordinate (s)",
        x_axis=time_bins,
        save=True,
        use_x_log=True,
        use_y_log=True
        )

def process_flux(transform):
    flux_scatter_data,raw_data= t.create_flux_scatter(modelFilePath, model_type, model, deltat=deltat, transform=transform,use_cache=True)
    time_bins = []
    for point in range(len(raw_data)):
        time_bins.append((point*step)+0.5)
    t.create_default_flux_plot(flux_scatter_data, "{model} Flux {transform}".format(model=model_type,transform=transform))
    t.create_regular_plot(
        plot_data=raw_data,
        x_axis=time_bins,
        axes_titles=[r'$\nu_x$', r'$\bar{\nu_e}$', r'$\nu_e$'],
        plot_title=f'{model_type} Truth Flux {transform}',
        ylab="Total Integrated Flux flavor/cm^2",
        xlab="Right Time in Coordinate (s)",
        use_x_log=True)

transforms_to_analyze = ['AdiabaticMSW_NMO','AdiabaticMSW_IMO']

for d in handlers.supported_detectors:
    for trans in transforms_to_analyze: process_detector(d,trans)

for trans in transforms_to_analyze: process_flux(trans)

# if multithreading==True:
#     for detector in profiles.keys():
#         proc = Process(target=process_detector, args=[detector])
#         proc.start()
#         proc.join()
# else:
#     for detector in detectors.keys():
#         process_detector(detector)