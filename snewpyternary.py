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

def create_detector_event_scatter(
        modelFilePath,
        model_type,
        detector,
        model,
        data_calc,
        deltat=1*u.s,
        d=10,
        transformation="NoTransformation",
        smearing=True,
        weighting="weighted"
        ):
    '''
    Creates normalized scatter data for use in a ternary diagram. Using SNOwGLoBES,
    a "truth flux" is created for a specified number of time bins. That is,
    in each time bin generated, a progenitor flux vs energy diagram is created.
    Then, this data is simmulated as though a detector would see the fluxes
    on Earth. Each time bin in the simulated output (now event rates vs energy)
    is integrated for the specified detector channel. From these detector channels,
    you must supply the data calculations since the data across detectors is not 
    uniform. To help with this, this method returns a dictionary of summed channel
    data to pick from.
    and then normalized against the other channels to create a ternary scatter point.
    

    Parameters
    ----------
    modelFilePath : str
        Path to the SNEWPY model. Must be an absolute path
    model_type : str
        SNEWPY name of the model
    detector : str
        SNOwGLoBES detector type. This detector must be available in $SNOWGLOBES
    model : snewpy.models
        The SNEWPY class abstraction of the model you're working with
    deltat : astropy.unit, optional
        Time bin width. The default is 1*u.s.
    d : int, optional
        simulated distance to progenitor. The default is 10.
    transformation : snewpy.transformation, optional
        Flavor transformation prescription. The default is "NoTransformation".
    smearing : bool, optional
        Use detector smearing matrix. The default is True.
    weighting : str "weighted" or "unweighted", optional
        Use detector data weights. The default is "weighted".
    data_calc : def
        Data calculation algorithm, since not all detector channel data is
        uniform. Passes a dictionary of available detector channels with their
        summed event rates in a given time bin

    Returns
    -------
    list
        The ternary scatter plot data

    '''
    
    if detector == 'all':
        raise("Cannot accept 'all' as detector type")
    snowglobes_out_name="snowglobes-output"
    snowglobes_dir = os.environ['SNOWGLOBES']
    
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
    print("\nNote that the data columns are as follows: a0, a1, a2, a2, ...")
    showing_columns = True
    for file in data_files: # which effectively goes through each time bin
        # we want to only have the ones that end in 'unweighted' for now
        if (weighting in file):
            # now fetch the time bin's energy chart
            e_bins = tables.get(file)
            
            if showing_columns:
                header = list(e_bins.get("header").split(" ")) # header in dictionary
                print(header)
                showing_columns = False
            
            # first we need to get the points by iterating through the points
            #tables.get(file).data[]
            data = e_bins.get("data") # 2D array for data
            
            dict_data = {}
            # build dictionary of available channels
            for i in range(len(header)):
                dict_data[header[i]]=np.sum(data[i])
            results = data_calc(dict_data)
            a= results[0]
            b=results[1]
            c=results[2]
            
            total = a+b+c
            plotting_data.append((100*a/total,100*b/total,100*c/total))
            
    # now retrieve header files
    #header_info = tables.get(data_files[1]).get("header").split(" ")
            
    return plotting_data

def create_default_detector_plot(plot_data,axes_titles,plot_title,save=True):
    '''
    From ternary detetector event scatter plot data, a ternary plot is created

    Parameters
    ----------
    plot_data : list of 3-tuples
        Scatter data produced by create_detector_event_scatter
    plot_title : str
        Plot title
    axes_titles : list
        List of axes titles in b,r,l order
    save : bool, optional
        Save the output file to ./plots/. The default is True.

    Returns
    -------
    figure, tax
        The figure and tax variables created by ternary.figure()

    '''
    figure, tax = ternary.figure(scale=100)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=10)
    tax.set_title(plot_title)
    # data is organized in top, right, left

    ### TODO: make sure that data_files[1] actually points to something that can get the header
    tax.bottom_axis_label(axes_titles[0])
    tax.right_axis_label(axes_titles[1])
    tax.left_axis_label(axes_titles[2])

    tax.scatter(points=plot_data, color="red")
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes

    tax.show()
    if save:
        tax.savefig('./plots/' + plot_title)
    return figure, tax
        
def create_flux_scatter(modelFilePath,
                        modeltype,
                        model,
                        deltat=1*u.s,
                        d=10,
                        transform="NoTransformation",
                        smearing=False):
    '''
    Similar to create_detector_event_scatter, although here we are just plotting
    the truth fluxes--the fluxes at the progenitor. A time series is created
    by SNOwGLoBES to make the flux vs energy graphs in time bins, then these
    fluxes are integrated along the time bins to create the normalized
    scatter plot data

    Parameters
    ----------
    modelFilePath : star
        Path to the SNEWPY model. Must be an absolute path.
    modeltype : str
        SNEWPY name of the model.
    model : snewpy.models
        The SNEWPY class abstraction of the model you’re working with.
    deltat : astropy.Quantity(), optional
        Time bin width. The default is 1*u.s.. The default is 1*u.s.
    d : int, optional
        Simulated distance to progenitor. The default is 10. The default is 10.
    transform : snewpy.transformation, optional
        Flavor transformation prescription. The default is “NoTransformation”. The default is "NoTransformation".
    smearing : bool, optional
        Use detector smearing matrix. The default is True. The default is False.

    Returns
    -------
    plotting_data : list
        The ternary scatter plot data.

    '''
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
        a=math.log(NuX+aNuX) if use_log else NuX
        b=math.log(aNuE) if use_log else aNuE
        c=math.log(NuE) if use_log else NuE
        total = a+b+c
        plotting_data.append((scale*a/total,scale*b/total,scale*c/total))
    return plotting_data

def create_default_flux_plot(plotting_data,plot_title,save=True):
    '''
    Creates a ternary plot from the truth flux data

    Parameters
    ----------
    plotting_data : list
        List of 3-tuples--the scatter points.
    plot_title : str
        The name of the plot.
    save : bool, optional
        Save the output file to ./plots/. The default is True.

    Returns
    -------
    figure, tax
        The figure and tax variables created by ternary.figure()

    '''
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
    return figure, tax