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
import sys
sys.path.insert(0,'./SURF2020')
from SURF2020.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

#simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
step = 0.04 #0.01 for best results
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

def process_detector(detector, transform):
# detector='scint20kt'
# transform='AdiabaticMSW_NMO'
    plot_data, raw_data, l_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                                detector,
                                                model,
                                                deltat=deltat,
                                                transformation=transform,
                                                data_calc=profiles[detector]['handler'],
                                                use_cache=True
                                                )
    # also create heatmap using Rishi's code
    heatmap_dict = generate_heatmap_dict(raw_data, plot_data)
    figure, tax = t.create_default_detector_plot(plot_data,
                                                  profiles[detector]['axes'](),
                                                  f'{model_type} {detector} {transform} Ternary',
                                                  heatmap=heatmap_dict,
                                                  show=False,
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
    return figure, tax, plot_data, raw_data, heatmap_dict

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

def remap_dict(dictionary,newval):
    # remaps a dictionary's 1 value to a different value
    new_dict = {}
    for k in dictionary.keys():
        if dictionary[k] == 1:
            new_dict[k] = newval
        else:
            new_dict[k] = 0
    return new_dict

# for d in handlers.supported_detectors:
for d in handlers.supported_detectors:
    nmo_figure, nmo_tax, nmo_plot_data, nmo_raw_data, nmo_hp = process_detector(d,transforms_to_analyze[0])
    imo_figure, imo_tax, imo_plot_data, imo_raw_data, imo_hp = process_detector(d,transforms_to_analyze[1])
    
    # just gonna have to create a new graph for combining things
    figure, tax = ternary.figure(scale=100)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="black", multiple=10)
    title = f'{model_type} {d} NMO vs IMO'
    tax.set_title(title)
    # data is organized in top, right, left
    
    axes_titles = profiles[d]['axes']()
    ### TODO: make sure that data_files[1] actually points to something that can get the header
    tax.bottom_axis_label(axes_titles[0])
    tax.right_axis_label(axes_titles[1])
    tax.left_axis_label(axes_titles[2])
    
    d1 = nmo_hp.copy()
    d2 = imo_hp.copy()
    newval = 0.6
    d2 = remap_dict(d2,newval)
    tax.heatmap(consolidate_heatmap_data(d1,d2),cbarlabel='In 68% CL')
    
    tax.scatter(points=nmo_plot_data, color='red')
    tax.scatter(points=imo_plot_data, color='blue')
    
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()
    # tax.get_axes().axis('off') # disables regular matlab plot axes
    tax.show()
    tax.savefig(f'./plots/{title}.png')

# for trans in transforms_to_analyze: process_flux(trans)

# if multithreading==True:
#     for detector in profiles.keys():
#         proc = Process(target=process_detector, args=[detector])
#         proc.start()
#         proc.join()
# else:
#     for detector in detectors.keys():
#         process_detector(detector)