#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 22:37:52 2022
Script used for plotting individual time slice data
@author: phyics
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.models import Nakazato_2013
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
import os
import ternary
import math

import io
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
import snewpyternary as t
import data_handlers as handlers

#simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
deltat=1*u.s
detector = 'wc100kt30prct'
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'
model
print(f'Timeframe of model is from {model.time[0]} to {model.time[len(model.time)-1]}')

transform = 'NoTransformation'
# create scintillator detector analysis
plot_data, raw_data, l_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                            detector,
                                            model,
                                            deltat=deltat,
                                            data_calc=handlers.h_wc100kt30prct,smearing=smearing)

# now iterate over select time slices
selected_bins = [0,5,10,15,19] #np.arange(0,19,step=1)#np.arange(0,model.time[-1],step=1)
for bin_no in selected_bins:
    nue_plus_es_energy = np.add(l_data[bin_no]['nue_O16'],l_data[bin_no]['e'])
    plt.plot(l_data[bin_no]['Energy'],nue_plus_es_energy,label="NuE+ES")
    plt.plot(l_data[bin_no]['Energy'],l_data[bin_no]['ibd'],label="ibd")
    plt.plot(l_data[bin_no]['Energy'],l_data[bin_no]['nc'],label="nc")
    
    # now also calculate event integrals
    ibd_integral = np.sum(l_data[bin_no]['ibd'])
    nue_plus_es_integral = np.sum(nue_plus_es_energy)
    nc_integral = np.sum(l_data[bin_no]['nc'])
    
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Event Rate')
    time_actual = (bin_no+1)-0.5 # in seconds, taken from the left side of the time bin
    title = f'{model_type} {detector} t={time_actual} (left in s)'
    plt.title(title)
    caption = f"Total Event Rate: {ibd_integral+nue_plus_es_integral+nc_integral}\nBin No: {bin_no}"
    y_max = float(list(plt.gca().get_ylim())[1])
    plt.text(0,-1*0.32*y_max,caption,ha='left')
    plt.legend()
    plt.savefig(f'./plots/time_slices/{title}.png',bbox_inches='tight')
    plt.show()


# now we take a time slice and plot its energy
figure, tax = t.create_default_detector_plot(plot_data,
                                              ['ibd','nue+es','nc'],
                                              '{model} {detector} {transform}'.format(model=model_type,detector=detector,transform=transform),
                                              save=True)
t.create_regular_plot(raw_data, ['ibd','nue+es','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector=detector,transform=transform),
                      ylab="Event Counts",save=True)