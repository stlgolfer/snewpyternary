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
import snewpyternary as t
import data_handlers as handlers
from matplotlib.animation import PillowWriter
import PIL

import sys
sys.path.insert(0,'./SURF2020')
from SURF2020.ternary_helpers import generate_heatmap_dict

#simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
step_size = 0.04
deltat=step_size*u.s
detector = 'wc100kt30prct'
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'
model
transform = "AdiabaticMSW_NMO"
print(f'Timeframe of model is from {model.time[0]} to {model.time[len(model.time)-1]}')

profiles = handlers.build_detector_profiles()
# create scintillator detector analysis
plot_data, raw_data, l_data = t.create_detector_event_scatter(
    modelFilePath,model_type,
    detector,
    model,
    deltat=deltat,
    transformation=transform,
    data_calc=profiles[detector]['handler'],
    use_cache=True
    )
heatmap_dict = generate_heatmap_dict(raw_data, plot_data)
figure, tax = t.create_default_detector_plot(plot_data,
                                              profiles[detector]['axes'](),
                                              f'{model_type} {detector} {transform} Ternary dt={str(deltat)}',
                                              show=True,
                                              heatmap=heatmap_dict,
                                              save=True)

# now iterate over select time slices
no_total_bins = int((model.time[-1].value-model.time[0].value)/step_size)
no_slices = 5
selected_bins = [] #np.arange(0,19,step=1)#np.arange(0,model.time[-1],step=1)
for x in range(no_slices):
    selected_bins.append(int(x*no_total_bins/no_slices))

# plot an bin with region as a ternary diagram
single_point_plots = []
# for bin_no in [0,5]:
bin_no=0
singular_normalized = plot_data[bin_no]
fig, ta = t.create_default_detector_plot(
    [singular_normalized],
    profiles[detector]['axes'](),
    f'Singular Bin {bin_no} Ternary dt={str(deltat)}',
    heatmap=generate_heatmap_dict([raw_data[bin_no]], [singular_normalized]),show=False,save=False
    )
fig.canvas.draw()
single_point_plots.append(plt.imshow(PIL.Image.frombytes('RGB',fig.canvas.get_width_height(),fig.canvas.tostring_rgb())))
# ta.savefig('test.png')

# here we will also create an animation of what just happened
ani = mpl.animation.ArtistAnimation(plt.figure(),single_point_plots)
ani.save('./plots/animation.gif', writer='pillow')

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
    time_actual = step_size*(bin_no+1)-0.5 # in seconds, taken from the left side of the time bin
    title = f'{model_type} {detector} t={time_actual} dt={deltat} (left in s)'
    plt.title(title)
    caption = f"Total Event Rate: {ibd_integral+nue_plus_es_integral+nc_integral}\nBin No: {bin_no}"
    y_max = float(list(plt.gca().get_ylim())[1])
    plt.text(0,-1*0.32*y_max,caption,ha='left')
    plt.legend()
    plt.savefig(f'./plots/time_slices/{title}.png',bbox_inches='tight')
    plt.show()

# now we plot as normal
figure, tax = t.create_default_detector_plot(plot_data,
                                              ['ibd','nue+es','nc'],
                                              '{model} {detector} {transform}'.format(model=model_type,detector=detector,transform=transform),
                                              save=True)
t.create_regular_plot(raw_data, ['ibd','nue+es','nc'],
                      '{model} {detector} {transform}'.format(model=model_type,detector=detector,transform=transform),
                      ylab="Event Counts",use_x_log=True,save=True)