# -*- coding: utf-8 -*-
"""
COPY OF ORIGINAL SNEWPY CODE. USED FOR ALTERNATIVE SNEWPY IMPLEMENTATIONS

The ``snewpy.snowglobes`` module contains functions for interacting with SNOwGLoBES.

`SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_ can estimate detected
event rates from a given input supernova neutrino flux. It supports many
different neutrino detectors, detector materials and interaction channels.
There are three basic steps to using SNOwGLoBES from SNEWPY:

* **Generating input files for SNOwGLoBES:**
    There are two ways to do this, either generate a time series or a fluence file. This is done taking as input the supernova simulation model.
    The first will evaluate the neutrino flux at each time step, the latter will compute the integrated neutrino flux (fluence) in the time bin.
    The result is a compressed .tar file containing all individual input files.
* **Running SNOwGLoBES:**
    This step convolves the fluence generated in the previous step with the cross-sections for the interaction channels happening in various detectors supported by SNOwGLoBES.
    It takes into account the effective mass of the detector as well as a smearing matrix describing the energy-dependent detection efficiency.
    The output gives the number of events detected as a function of energy for each interaction channel, integrated in a given time window (or time bin), or in a snapshot in time.
* **Collating SNOwGLoBES outputs:**
    This step puts together all the interaction channels and time bins evaluated by SNOwGLoBES in a single file (for each detector and for each time bin).
    The output tables allow to build the detected neutrino energy spectrum and neutrino time distribution, for each reaction channel or the sum of them.
"""

import io
import logging
import os
import re
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from tqdm.auto import tqdm
from warnings import warn
import math

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy
#from snewpy.snowglobes_interface import SNOwGLoBES, SimpleRate

logger = logging.getLogger(__name__)

def w_generate_time_series(model_path,
                         model_type,
                         transformation_type,
                         d,
                         output_filename=None,
                         ntbins=30,
                         deltat=None,
                         snmodel_dict={},
                         log_bins=False):
    """Generate time series files in SNOwGLoBES format.

    This version will subsample the times in a supernova model, produce energy
    tables expected by SNOwGLoBES, and compress the output into a tarfile.

    Parameters
    ----------
    model_path : str
        Input file containing neutrino flux information from supernova model.
    model_type : str
        Format of input file. Matches the name of the corresponding class in :py:mod:`snewpy.models`.
    transformation_type : str
        Name of flavor transformation. See snewpy.flavor_transformation documentation for possible values.
    d : int or float
        Distance to supernova in kpc.
    output_filename : str or None
        Name of output file. If ``None``, will be based on input file name.
    ntbins : int
        Number of time slices. Will be ignored if ``deltat`` is also given.
    deltat : astropy.Quantity or None
        Length of time slices.
    snmodel_dict : dict
        Additional arguments for setting up the supernova model. See documentation of relevant ``SupernovaModel`` subclass for available options. (Optional)
    log_bins : bool
        Use log bins or not

    Returns
    -------
    str
        Path of compressed .tar file with neutrino flux data.
    """
    model_class = getattr(snewpy.models.ccsn, model_type)

    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    model_dir, model_file = os.path.split(os.path.abspath(model_path))
    snmodel = model_class(model_path, **snmodel_dict)

    # Subsample the model time. Default to 30 time slices.
    tmin = snmodel.get_time()[0]
    tmax = snmodel.get_time()[-1]
    if deltat is not None:
        dt = deltat
        ntbins = int((tmax-tmin)/dt)
    else:
        dt = (tmax - tmin) / (ntbins+1)

    tedges = np.arange(tmin/u.s, tmax/u.s, dt/u.s)*u.s
    times = 0.5*(tedges[1:] + tedges[:-1])
    # print(list(times))
    
    # now process log data
    if (log_bins==True):
        needed_offset = -1*u.s
        # need to offset the times so no negatives
        if tedges[0] < 0:
            needed_offset = tedges[0] + 0.0001*u.s
            tedges = tedges + tedges[0] + 0.0001*u.s #shift so it's very close to 0
            tmin = 0.0001*u.s
            tmax+=tmin
        log_edges = np.asarray([])
        
        tstep = math.log10(abs(tmax/tmin))/len(times)
        for i in range(0,len(times)):
            t = (tmin/u.s)*(10**(i*tstep))
            # print(t)
            log_edges = np.append(log_edges,t)
        # log_edges = np.logspace(math.log((tmin/u.s).value),math.log((tmax/u.s).value),len(times))
        # print((tmax/u.s).value)
        log_edges = log_edges*u.s
        if needed_offset>-1*u.s:
            #keep in mind needed_offset will still be a negative number here 
            log_edges = log_edges + needed_offset
        times = 0.5*(log_edges[1:] + log_edges[:-1])
        print(f'Proceeding with {len(times)} bin(s)')
        # plt.figure()
        # plt.scatter([no for no in range(len(times))], times)
        # plt.show()
        # print(len(times))

    # Generate output.
    if output_filename is not None:
        tfname = output_filename + 'kpc.tar.bz2'
    else:
        model_file_root, _ = os.path.splitext(model_file)  # strip extension (if present)
        tfname = model_file_root + '.' + transformation_type + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin, tmax, ntbins, d) + 'kpc.tar.bz2'

    with tarfile.open(os.path.join(model_dir, tfname), 'w:bz2') as tf:
        #creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str, transformation_type)).encode('ascii')
        tf.addfile(tarfile.TarInfo(name='parameterinfo'), io.BytesIO(output))

        MeV = 1.60218e-6 * u.erg
        energy = np.linspace(0, 100, 501) * MeV  # 1MeV

        # Loop over sampled times.
        for i, t in enumerate(times):
            osc_spectra = snmodel.get_transformed_spectra(t, energy, flavor_transformation)

            osc_fluence = {}
            table = []

            table.append('# TBinMid={:g}sec TBinWidth={:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2'.format(t, dt))
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 0.2 * MeV / (4.*np.pi*(d*1000*3.086e+18)**2)

                s = '{:17.8E}'.format(E/(1e3 * MeV))
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            model_file_root, _ = os.path.splitext(model_file)
            filename = model_file_root + '.tbin{:01d}.'.format(i+1) + transformation_type + \
                '.{:.3f},{:.3f},{:01d}-{:.1f}kpc{}'.format(tmin/u.s, tmax/u.s, ntbins, d, extension)

            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))

    return os.path.join(model_dir, tfname)
