#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 00:15:52 2022
here we will need to open the computed fluxes and then integrate over the
fluence over all energy.
We will use snowglobes to generate the fluence files, then we will integrate
over all energies to get total flux for that time bin for each flavor. then
plot as before
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

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

# simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
#ntbins =  # number of time bins
deltat = (1*u.s) # time bin size, details found in SNEWPY article
d = 10 # in pc, distance to SN

snowglobes_dir = os.environ['SNOWGLOBES']
tball_suffix = 'kpc.tar.bz2'
#print(os.environ['SNOWGLOBES'])
smearing = True
model

'''
# we have the model loaded, so let's try to generate the fluence files using SNEWPY
snowglobes.generate_fluence(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=10,
                            output_filename=snowglobes_out_name
                            )
'''

tarball_path = snowglobes.generate_time_series(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=d,
                            deltat=deltat
                            )
fluence_data = []
labels_data = None
with TemporaryDirectory(prefix='snowglobes') as tempdir:
    with tarfile.open(tarball_path) as tar:
        tar.extractall(tempdir)

    flux_files = list(Path(tempdir).glob('*.dat'))
    for flux in flux_files:
        fluence_data.append(np.loadtxt(str(flux),unpack=True))
        print('NEW FLUENCE\n================\n'+str(flux)+'\n')
        f=open(str(flux), "r")
        print(f.readline()) # prints .dat data file header info
        labels = f.readline()
        print(labels) # prints .dat file info
        
        # generate labels_data if not already; uniformity in files means we
        # only need to do this once
        if (labels_data == None):
            labels_data=labels.split("	")
        f.close()
# now that we have all the fluences loaded per time bin, we now need to
# integrate the fluence to get the total flux
# data comes out backwards, so first need to flip it
fluence_data.reverse() # now they're in the correct time sequence

plotting_data = []
for time_bin in fluence_data:
    NuE = np.sum(time_bin[1])
    NuX = np.sum(time_bin[2])+np.sum(time_bin[3])
    aNuE = np.sum(time_bin[4])
    aNuX = np.sum(time_bin[5])+np.sum(time_bin[6])
    a=NuX
    b=aNuE
    c=NuE
    total = a+b+c
    plotting_data.append((100*a/total,100*b/total,100*c/total))
    
plot_title = "Nakazato_2013 Fluxes"
figure, tax = ternary.figure(scale=100)
tax.boundary(linewidth=2.0)
tax.gridlines(color="blue", multiple=20)
tax.set_title(plot_title)
# data is organized in top, right, left

### TODO: make sure that data_files[1] actually points to something that can get the header
tax.bottom_axis_label(r'$\nu_x$')
tax.right_axis_label(r'$\bar{\nu_e}$')
tax.left_axis_label(r'$\nu_e$')

tax.scatter(points=plotting_data, color="red")
tax.ticks(axis='lbr', linewidth=1, multiple=20)
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off') # disables regular matlab plot axes

tax.show()
tax.savefig('./plots/' + plot_title)