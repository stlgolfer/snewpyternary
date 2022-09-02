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
import multiprocessing as mp
import caching as cache
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data
from model_wrappers import snewpy_models, sn_model_default_time_step

#simulation details
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'

# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
transform_list = list(flavor_transformation_dict.keys())
transforms_to_analyze = ['NoTransformation'] #['AdiabaticMSW_NMO','AdiabaticMSW_IMO','NoTransformation']
profiles = handlers.build_detector_profiles()

# param parser
if len(sys.argv) <= 1:
    raise RuntimeError("Usage: python meta_analysis.py --show=<True | False>")

show_charts: bool = True if sys.argv[1].split('--show=')[1] == "True" else False

def process_detector(config: t.MetaAnalysisConfig, detector: str) -> None:
    plot_data, raw_data, l_data = t.create_detector_event_scatter(
        config.model_file_path,
        config.model_type,
        detector,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transformation=config.transformation,
        data_calc=profiles[detector]['handler'],
        use_cache=True,
        log_bins=True
    )
    # also create heatmap using Rishi's code
    # heatmap_dict = generate_heatmap_dict(raw_data, plot_data)
    figure, tax = t.create_default_detector_plot(plot_data,
                                                  profiles[detector]['axes'](),
                                                  f'{config.model_type} {detector} {config.transformation} Ternary',
                                                  show=show_charts,
                                                  save=True)
    return plot_data, raw_data

def process_flux(config: t.MetaAnalysisConfig) -> None:
    flux_scatter_data,raw_data= t.create_flux_scatter(
        config.model_file_path,
        config.model_type,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transform=config.transformation,
        use_cache=True,
        log_bins=True
    )
    t.create_default_flux_plot(
        flux_scatter_data,
        "{model} Flux {transform}".format(model=config.model_type,transform=config.transformation),
        show=show_charts
        )
    
    t.create_regular_plot(
        plot_data=raw_data,
        axes_titles=[r'$\nu_x$', r'$\bar{\nu_e}$', r'$\nu_e$'],
        plot_title=f'{config.model_type} Truth Flux {config.transformation}',
        ylab="Total Integrated Flux flavor/cm^2",
        xlab="Right Time in Coordinate (s)",
        show=show_charts,
        use_x_log=False,save=True)

def remap_dict(dictionary,newval):
    # remaps a dictionary's 1 value to a different value
    new_dict = {}
    for k in dictionary.keys():
        if dictionary[k] == 1:
            new_dict[k] = newval
        else:
            new_dict[k] = 0
    return new_dict

def process_transformation(config: t.MetaAnalysisConfig):
    print(f'Now processing {config.model_type}')
    # first get the flux data
    process_flux(config)
    
    p_data, r_data = process_detector(config,'ar40kt')
    # need to convert data to an array
    all_plot_data = [list(key) for key in r_data]# going to take each detector and add them up
    
    for detector in ['wc100kt30prct','scint20kt']:
        p_data, r_data = process_detector(config,detector)
        all_plot_data = all_plot_data + np.asarray([list(key) for key in r_data])
    # need to figure out a way to sum all the detectors
    # now renormalize and convert all points back to tuples
    t.create_regular_plot(all_plot_data, handlers.same_axes(), 'Regular plot of *Detectors', ylab='Event rate',show=show_charts)
    
    normalized = []
    for point in all_plot_data:
        a=point[0]
        b=point[1]
        c=point[2]
        tot=a+b+c
        normalized.append((100*a/tot,100*b/tot,100*c/tot))
    # all_plot_data = [tuple(point[0]) for point in all_plot_data]
    t.create_regular_plot(normalized, handlers.same_axes(), 'Super normalized ternary points', 'Event Rate',show=show_charts)
    
    scale=100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale/10)
    title=f'{config.model_type} *Detectors {config.transformation} Ternary Logged Bins'
    tax.set_title(title)
    # data is organized in top, right, left

    tax.bottom_axis_label('nux')
    tax.right_axis_label('nuebar')
    tax.left_axis_label('nue')

    tax.scatter(points=normalized, color="red")
    tax.ticks(axis='lbr', linewidth=1, multiple=scale/10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes
    tax.savefig(f'./all_detector_plots/{title}')
    
    if show_charts == True:
        tax.show()

# process_transformation(t.MetaAnalysisConfig(snewpy_models['Bollig_2016'], 'NoTransformation'))

if __name__ == '__main__': # for multiprocessing
    for transformation in transforms_to_analyze:
        for model in snewpy_models.keys():
            proc = mp.Process(target=process_transformation, args=[t.MetaAnalysisConfig(snewpy_models[model], transformation)])
            proc.start()
            proc.join()

# for d in handlers.supported_detectors:
# for d in handlers.supported_detectors:
#     nmo_figure, nmo_tax, nmo_plot_data, nmo_raw_data, nmo_hp = process_detector(d,transforms_to_analyze[0])
#     imo_figure, imo_tax, imo_plot_data, imo_raw_data, imo_hp = process_detector(d,transforms_to_analyze[1])
    
#     # just gonna have to create a new graph for combining things
#     figure, tax = ternary.figure(scale=100)
#     tax.boundary(linewidth=2.0)
#     tax.gridlines(color="black", multiple=10)
#     title = f'{model_type} {d} NMO vs IMO'
#     tax.set_title(title)
#     # data is organized in top, right, left
    
#     axes_titles = profiles[d]['axes']()
#     ### TODO: make sure that data_files[1] actually points to something that can get the header
#     tax.bottom_axis_label(axes_titles[0])
#     tax.right_axis_label(axes_titles[1])
#     tax.left_axis_label(axes_titles[2])
#     d1 = nmo_hp.copy()
#     d2 = imo_hp.copy()
#     newval = 0.6
#     d2 = remap_dict(d2,newval)
#     tax.heatmap(consolidate_heatmap_data(d1,d2),cbarlabel='In 68% CL')
    
#     tax.scatter(points=nmo_plot_data, color='red')
#     tax.scatter(points=imo_plot_data, color='blue')
    
#     tax.ticks(axis='lbr', linewidth=1, multiple=10)
#     tax.clear_matplotlib_ticks()
#     # tax.get_axes().axis('off') # disables regular matlab plot axes
#     tax.show()
#     tax.savefig(f'./plots/{title}.png')

# for trans in transforms_to_analyze: process_flux(trans)

# if multithreading==True:
#     for detector in profiles.keys():
#         proc = Process(target=process_detector, args=[detector])
#         proc.start()
#         proc.join()
# else:
#     for detector in detectors.keys():
#         process_detector(detector)