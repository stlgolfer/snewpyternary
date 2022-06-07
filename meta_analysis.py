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
model

detector_list = ['scint20kt','ar40kt','wc100kt30prct']
# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
transform_list = list(flavor_transformation_dict.keys())

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

transform = 'NoTransformation'
# create scintillator detector analysis
plot_data, raw_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                            'scint20kt',
                                            model,
                                            data_calc=h_scint20kt)
figure, tax = t.create_default_detector_plot(plot_data,
                                              ['ibd','nue+es','nc'],
                                              '{model} {detector} {transform}'.format(model=model_type,detector='scint20kt',transform=transform),
                                              save=True)
t.create_regular_plot(raw_data, ['ibd','nue+es','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector='scint20kt',transform=transform),
                      ylab="Event Counts",save=True)


# create argon detector analysis
plot_data, raw_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                            'ar40kt',
                                            model,
                                            data_calc=h_ar40kt)
t.create_default_detector_plot(plot_data,
                                ['nue','nuebar','nc'],
                                '{model} {detector} {transform}'.format(model=model_type,detector='ar40kt',transform=transform))
t.create_regular_plot(raw_data, ['nue','nuebar','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector='ar40kt',transform=transform),
                      ylab="Event Counts",save=True)

# create water detector analysis
plot_data, raw_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                            'wc100kt30prct',
                                            model,
                                            data_calc=h_wc100kt30prct)
t.create_default_detector_plot(plot_data,
                                ['ibd','nue+es','nc'],
                                '{model} {detector} {transform}'.format(model=model_type,detector='wc100kt30prct',transform=transform))
t.create_regular_plot(raw_data,
                      ['ibd','nue+es','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector='wc100kt30prct',transform=transform),
                      ylab="Event Counts",save=True)


# for saving figures with captions
# fig = tax.get_figure()
# fig.text(x=0.5,y=0.001,s="Some subtext goes here")
# fig.savefig(fname='./plots/test.png')

# for d in detector_list:
#     # for each detector config, simulate all types of transformations
#     for transform in transform_list:
#         plot_title = "Nakazato_2013 Events " + weighting + " "  + d + " " + transform
#         plot_data = t.create_detector_event_scatter(modelFilePath,model_type, d, model,transformation=transform)
#         # then we can generate a ternary plot for this
#         t.create_default_detector_plot(plot_data,plot_title)

# create flux plot as well
# for transform in transform_list:
#     flux_scatter_data = t.create_flux_scatter(modelFilePath, model_type, model,transform=transform)
#     t.create_default_flux_plot(flux_scatter_data, "{model} Flux {transform}".format(model=model_type,transform=transform))
flux_scatter_data,raw_data = t.create_flux_scatter(modelFilePath, model_type, model,transform=transform)
t.create_default_flux_plot(flux_scatter_data, "{model} Flux {transform}".format(model=model_type,transform=transform))
fig = t.create_regular_plot(raw_data, ['NuX', 'aNuE', 'NuE'], f'{model_type} Truth Flux {transform}', ylab="Total Integrated Flux flavor/cm^2")
# fig.savefig(f'./plots/{model_type} Truth Flux {transform}')
# fig.show()