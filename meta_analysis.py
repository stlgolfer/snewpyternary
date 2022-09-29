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
from ternary import TernaryAxesSubplot

import snewpyternary as t
import os
import ternary
import math
from snewpy.flavor_transformation import *
import data_handlers as handlers
import multiprocessing as mp
import caching as cache
import sys

import snowglobes_wrapper

sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data
from model_wrappers import snewpy_models, sn_model_default_time_step
import click

#simulation details
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'

# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
complete_transform_list = list(flavor_transformation_dict.keys())
transforms_to_analyze = complete_transform_list # ['NoTransformation'] #['AdiabaticMSW_NMO','AdiabaticMSW_IMO','NoTransformation']
profiles = handlers.build_detector_profiles()

# global params
show_charts: bool = True
use_log: bool = True
use_cache: bool = True

_colors = ['RED', 'GREEN', 'BLUE']

def process_detector(config: t.MetaAnalysisConfig, set_no: int, detector: str) -> None:
    plot_data, raw_data, l_data = t.create_detector_event_scatter(
        config.model_file_paths[set_no],
        config.model_type,
        detector,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transformation=config.transformation,
        data_calc=profiles[detector]['handler'],
        use_cache=use_cache,
        log_bins=use_log
    )
    # also create heatmap using Rishi's code
    # heatmap_dict = generate_heatmap_dict(raw_data, plot_data)
    figure, tax = t.create_default_detector_plot(plot_data,
                                                  profiles[detector]['axes'](),
                                                  f'{config.model_type} {detector} {config.transformation} Ternary',
                                                  show=show_charts,
                                                  save=True)
    return plot_data, raw_data

def process_flux(config: t.MetaAnalysisConfig, set_no: int) -> None:
    flux_scatter_data,raw_data = t.create_flux_scatter(
        config.model_file_paths[set_no],
        config.model_type,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transform=config.transformation,
        use_cache=use_cache,
        log_bins=use_log
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

def aggregate_detector(config: t.MetaAnalysisConfig, number: int, colorid: int, tax: TernaryAxesSubplot) -> None:
    process_flux(config, number)

    # print out information of the set
    print(config.model(config.model_file_paths[number]))

    p_data, r_data = process_detector(config, number, 'ar40kt')
    # need to convert data to an array
    all_plot_data = [list(key) for key in r_data]  # going to take each detector and add them up

    for detector in ['wc100kt30prct', 'scint20kt']:
        p_data, r_data = process_detector(config, number, detector)
        all_plot_data = all_plot_data + np.asarray([list(key) for key in r_data])

    # now get the time bins
    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[number],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=use_log
    )

    t.create_regular_plot(all_plot_data,
                          handlers.same_axes(),
                          f'*Detectors {config.model_type} {config.transformation} {_colors[number]} {config.model_file_paths[number].split("/")[-1]}.png',
                          x_axis=time_bins_x_axis,
                          ylab='Event rate',
                          show=show_charts
                          )

    # now renormalize and convert all points back to tuples
    normalized = []
    for point in all_plot_data:
        a = point[0]
        b = point[1]
        c = point[2]
        tot = a + b + c
        normalized.append((100 * a / tot, 100 * b / tot, 100 * c / tot))
    # all_plot_data = [tuple(point[0]) for point in all_plot_data]
    t.create_regular_plot(normalized, handlers.same_axes(), f'{config.model_type} Super Normalized Ternary Points', 'Event Rate',
                          show=show_charts)

    # going to try dynamically sized points between lines?
    widths = np.linspace(0.01, 1, num=len(normalized))
    for p in range(len(normalized) - 1):
        if (p + 1 >= len(normalized)):
            break
        tax.line(normalized[p], normalized[p + 1], color=(widths[p] if colorid == 0 else 0, widths[p] if colorid == 1 else 0, widths[p] if colorid == 2 else 0, 1), linestyle=':', linewidth=3)


def process_transformation(config: t.MetaAnalysisConfig):
    print(f'Now processing {config.model_type}')
    # first get the flux data

    
    scale=100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale/10)
    title=f'{config.model_type} *Detectors {config.transformation} Logged Bins'
    tax.set_title(title)
    # data is organized in top, right, left

    tax.bottom_axis_label('nux')
    tax.right_axis_label('nuebar')
    tax.left_axis_label('nue')

    # heatmap line goes here
    #timemap = {}
    #colorspace = [c for c in range(len(all_plot_data))]
    #for j, n in enumerate(normalized):
        #timemap[n] = float(j+10)
        # then vary each coordinate
        #p_i = n # initial tuple coordinate in ternary space
        #for r in range(heatmap_radius):
            #p_pos = (int(p_i[0]) + r,int(p_i[1]) + r, int(p_i[2]) +r)
            #p_neg = (int(p_i[0]) - r,int(p_i[1]) - r, int(p_i[2]) -r)
            # add to timemap
            #timemap[p_pos] = j
            #timemap[p_neg] = j
            
    #if show_time_heatmap == True:
    #    tax.heatmap(timemap)

    # need to create different colors
    colorid: int = 0
    f = open(f'./all_detector_plots/{title}.txt', 'w')
    for set_number in config.set_numbers:
        # write model metadata as well
        f.write(f'{_colors[colorid]}\n==========\n')
        f.write(repr(config.model(config.model_file_paths[set_number])))
        f.write('\n\n')
        aggregate_detector(config,set_number, colorid, tax)
        colorid+=1
    f.close()

    #tax.scatter(points=normalized)
    
    tax.ticks(axis='lbr', linewidth=1, multiple=scale/10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes
    tax.savefig(f'./all_detector_plots/{title}')
    
    if show_charts == True:
        tax.show()

# process_transformation(t.MetaAnalysisConfig(snewpy_models['Bollig_2016'], 'NoTransformation'))
@click.command()
@click.argument('models',required=True,type=str,nargs=-1)
@click.option('--showc',default=False,type=bool,help='Whether to show generated plots or not. Will always save and cache')
@click.option('-p',required=False, multiple=True, type=str, default=['NoTransformation'], help='Prescriptions to use')
@click.option('--distance',default=10,type=int,help='The distance (in kPc) to the progenitor source')
@click.option('--uselog',default=True,type=bool, help='Use logarithmically spaced (and shifted) time bins')
@click.option('--setno', required=False, default=[0],type=int,multiple=True, help='Model set index. See model_wrappers.py')
@click.option('--cache', required=False, default=True, type=bool, help='If true, use cache')
def start(showc,models,distance,uselog,p, setno, cache):
    global show_charts
    show_charts = showc
    
    global d
    d = distance
    
    global use_log
    use_log = uselog

    global use_cache
    use_cache = cache

    # check set numbers
    if len(setno) > 3:
        raise ValueError("Can only superimpose a maximum of 3 sets onto one chart")

    # want to iterate by model, then prescription, then set

    for model in models:
        for prescription in (flavor_transformation_dict.keys() if model == "ALL" else p):
            # check to see if valid model set number

            for no in setno:
                if no >= len(snewpy_models[model].file_paths):
                    raise ValueError(f"Invalid model set id. Max is {len(snewpy_models[model].file_paths)-1}")

            # remove multithreading for now. run sequentially
            process_transformation(t.MetaAnalysisConfig(snewpy_models[model], setno, prescription))

            # proc = mp.Process(target=process_transformation, args=[])
            # proc.start()
            # proc.join()
            
if __name__ == '__main__': # for multiprocessing
    start()

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
