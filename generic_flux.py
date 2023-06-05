import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.models import Nakazato_2013, OConnor_2015
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
from ternary import TernaryAxesSubplot

import data_handlers
import snewpyternary as t
import os
import ternary
import math
from snewpy.flavor_transformation import *
import data_handlers as handlers
import multiprocessing as mp
import caching as cache
import sys
import pickle
from tqdm import tqdm
from scipy.integrate import simpson
import Rebinning
from scipy.stats import norm

import snowglobes_wrapper

sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data
from model_wrappers import snewpy_models, sn_model_default_time_step
import click

def start():
    print("Starting")
    # first make the config
    config = t.MetaAnalysisConfig(
        snewpy_models['Nakazato_2013'],
        0,
        'AdiabaticMSW_IMO',
        proxy_config=data_handlers.ConfigBestChannel()
    )

    flux_scatter_data, raw_data, labeled = t.create_flux_scatter(
        config.model_file_paths[0],
        config.model_type,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transform=config.transformation,
        use_cache=True,
        log_bins=True
    )

    raw_data_restriped = np.array(list(zip(*raw_data)))

    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[0],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=True
    )

    fx_plot, fx_axes = plt.subplots(1, 1, figsize=(8, 8))
    truth_flux_labels = config.proxyconfig.flux_axes()
    fx_axes.plot(time_bins_x_axis, raw_data_restriped[0], linestyle='None', marker='.', label=truth_flux_labels[0])
    # fx_axes.plot(time_bins_x_axis, raw_data_restriped[1], linestyle='None', marker='.', label=truth_flux_labels[1])
    # fx_axes.plot(time_bins_x_axis, raw_data_restriped[2], linestyle='None', marker='.', label=truth_flux_labels[2])

    # gaussian fit for just one flavor (flavor 0)
    mean, std = norm.fit(raw_data_restriped[0])
    print(f'mean: {mean}, std: {std}')
    fx_axes.plot(time_bins_x_axis, norm.pdf(time_bins_x_axis, mean, std), label='fit', color='red')

    fx_axes.set_xlabel('Time (s)')
    fx_axes.set_ylabel(r'$neutrinos/cm^2$')
    fx_title = f'Truth Flux for \n{config.model_file_paths[0].split("/")[-1]}'
    fx_axes.set_title(fx_title)
    fx_axes.set_xscale('log')
    fx_axes.legend()
    fx_plot.savefig(f'./generic_flux/{t.clean_newline(fx_title)}.png')


if __name__ == '__main__':
    start()
