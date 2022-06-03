# goal of this preliminary chart will be to generate fluences given a model
# then create a ternary model of fluences and then also the event rates
# first, let's import the snewpy stuff
# let's start with the Nakazato model
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

'''
Creates normalized ternary scatter plot data from a snowglobes simulation
'''
def create_detector_event_scatter(
        modelFilePath,
        model_type,
        detector,
        model,
        deltat=1*u.s,
        d=10,
        transformation="NoTransformation",
        smearing=False,
        weighting="unweighted"
        ):
    
    if detector == 'all':
        raise("Cannot accept 'all' as detector type")
    # simulation details
    #modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
    #modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
    #model = Nakazato_2013(modelFilePath)
    #deltat = (1*u.s) # time bin size, details found in SNEWPY article
    #d = 10 # in pc, distance to SN
    #detector = "scint20kt" # refrain from using "all"
    snowglobes_out_name="snowglobes-output"
    snowglobes_dir = os.environ['SNOWGLOBES']
    #print(os.environ['SNOWGLOBES'])
    #smearing = True
    #model
    
    '''
    # we have the model loaded, so let's try to generate the fluence files using SNEWPY
    snowglobes.generate_fluence(model_path=modelFilePath,model_type="Nakazato_2013",
                                transformation_type="NoTransformation",
                                d=10,
                                output_filename=snowglobes_out_name
                                )
    '''
    
    tball_complete = snowglobes.generate_time_series(
        model_path=modelFilePath,
        model_type=model_type,
        transformation_type=transformation,
        d=d,
        output_filename=snowglobes_out_name,
        deltat=deltat
        )
    print("NEW SNOwGLoBES Simulation\n=============================")
    print("Using detector schema: " + detector)
    print("Using transform: " + transformation)
    print("=============================")
    
    # also need to run the simulation so that we get the correct output
    snowglobes.simulate(SNOwGLoBESdir=snowglobes_dir,
                        tarball_path=tball_complete,
                        detector_effects=smearing,
                        detector_input=detector
                        )
    
    '''
    this method is too collated--as a first try we just need the event calculations
    '''
    tables = snowglobes.collate(tarball_path=tball_complete,
                       smearing=smearing,skip_plots=True,SNOwGLoBESdir=snowglobes_dir)
    
    data_files = list(tables.keys())
    plotting_data = []
    
    '''
    so what we have to do is go through each time bin, sum each particle's event
    count, then append to plotting_data. Then, we will have a list of time-evolved
    count
    '''
    
    for file in data_files: # which effectively goes through each time bin
        # we want to only have the ones that end in 'unweighted' for now
        if (weighting in file):
            # now fetch the time bin's energy chart
            e_bins = tables.get(file)
            
            # first we need to get the points by iterating through the points
            header = list(e_bins.get("header").split(" ")) # header in dictionary
            print(header)
            #tables.get(file).data[]
            data = e_bins.get("data") # 2D array for data
            # now go through each available particles
            a=np.sum(data[1])
            b=np.sum(data[2])
            c=np.sum(data[3])
            total = a+b+c
            plotting_data.append((100*a/total,100*b/total,100*c/total))
            
    # now retrieve header files
    header_info = tables.get(data_files[1]).get("header").split(" ")
            
    return [plotting_data, [header_info[1],header_info[2],header_info[3]]]

'''
'''
def create_default_detector_plot(plot_data,plot_title,save=True):
    figure, tax = ternary.figure(scale=100)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=10)
    tax.set_title(plot_title)
    # data is organized in top, right, left

    ### TODO: make sure that data_files[1] actually points to something that can get the header
    tax.bottom_axis_label(plot_data[1][0])
    tax.right_axis_label(plot_data[1][1])
    tax.left_axis_label(plot_data[1][2])

    tax.scatter(points=plot_data[0], color="red")
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes

    tax.show()
    if save:
        tax.savefig('./plots/' + plot_title)
        
def create_flux_scatter(modelFilePath,
                        modeltype,
                        model,
                        deltat=1*u.s,
                        d=10,
                        transform="NoTransformation",
                        smearing=False):
    print("NEW SNOwGLoBES FLUENCE TIME SERIES GENERATION\n=============================")
    #print("Using detector schema: " + detector)
    print("Using transform: " + transform)
    print("=============================")
    tarball_path = snowglobes.generate_time_series(
        model_path=modelFilePath,
        model_type=modeltype,
        transformation_type=transform,
        d=d,
        deltat=deltat
        )
    fluence_data = []
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

    # now that we have all the fluences loaded per time bin, we now need to
    # integrate the fluence to get the total flux
    # data comes out backwards, so first need to flip it
    fluence_data.reverse() # now they're in the correct time sequence
    scale = 100
    use_log = False
    plotting_data = []
    for time_bin in fluence_data:
        NuE = np.sum(time_bin[1])
        NuX = np.sum(time_bin[2])+np.sum(time_bin[3])
        aNuE = np.sum(time_bin[4])
        aNuX = np.sum(time_bin[5])+np.sum(time_bin[6])
        a=math.log(NuX) if use_log else NuX
        b=math.log(aNuE) if use_log else aNuE
        c=math.log(NuE) if use_log else NuE
        total = a+b+c
        plotting_data.append((scale*a/total,scale*b/total,scale*c/total))
    return plotting_data

def create_default_flux_plot(plotting_data,plot_title,save=True):
    scale=100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale/10)
    tax.set_title(plot_title)
    # data is organized in top, right, left

    ### TODO: make sure that data_files[1] actually points to something that can get the header
    tax.bottom_axis_label(r'$\nu_x$')
    tax.right_axis_label(r'$\bar{\nu_e}$')
    tax.left_axis_label(r'$\nu_e$')

    tax.scatter(points=plotting_data, color="red")
    tax.ticks(axis='lbr', linewidth=1, multiple=scale/10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes

    tax.show()
    tax.savefig('./plots/' + plot_title)