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

def h_scint20kt(data):
    # must return a list of a, b, c
    ibd = np.sum(data['ibd'])
    nue_plus_es=np.sum(data['nue_C12'])+np.sum(data['nue_C13']+data['e'])
    nc = np.sum(data['nc'])
    return [ibd,nue_plus_es,nc]

transform = 'NoTransformation'
# create scintillator detector analysis
plot_data, raw_data, l_data = t.create_detector_event_scatter(modelFilePath,model_type,
                                            'scint20kt',
                                            model,
                                            data_calc=h_scint20kt)

# now iterate over select time slices
selected_bins = [0, 3, 10, 12, 18]
for bin_no in selected_bins:
    nue_plus_es_energy = np.add(np.add(l_data[bin_no]['nue_C12'],l_data[bin_no]['nue_C13']),l_data[bin_no]['e'])
    plt.plot(l_data[bin_no]['Energy'],nue_plus_es_energy,label="NuE+ES")
    plt.plot(l_data[bin_no]['Energy'],l_data[bin_no]['ibd'],label="ibd")
    plt.plot(l_data[bin_no]['Energy'],l_data[bin_no]['nc'],label="nc")
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Flux at Earth in flavor/cm^2')
    time_actual = (bin_no+1)-0.5 # in seconds, taken from the left side of the time bin
    title = f'{model_type} scint20kt t={time_actual} (left in s)'
    plt.title(title)
    plt.legend()
    plt.savefig(f'./plots/time_slices/{title}.png')
    plt.show()


# now we take a time slice and plot its energy
figure, tax = t.create_default_detector_plot(plot_data,
                                              ['ibd','nue+es','nc'],
                                              '{model} {detector} {transform}'.format(model=model_type,detector='scint20kt',transform=transform),
                                              save=True)
t.create_regular_plot(raw_data, ['ibd','nue+es','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector='scint20kt',transform=transform),
                      ylab="Event Counts",save=True)