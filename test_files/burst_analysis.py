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
import numpy as np

import snowglobes_wrapper

sys.path.insert(0,'../SURF2020fork')
from SURF2020fork.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data
from model_wrappers import snewpy_models, sn_model_default_time_step, SNEWPYModel
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

def process_detector(config: t.MetaAnalysisConfig, set_no: int, detector: str):
    plot_data, raw_data, l_data = t.create_detector_event_scatter(
        config.model_file_paths[set_no],
        config.model_type,
        detector,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transformation=config.transformation,
        data_calc=profiles[detector]['handler'],
        use_cache=True
    )
    return plot_data, raw_data

def aggregate_detector(config: t.MetaAnalysisConfig, number: int) -> None:
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
        log_bins=False,
        presn=False
    )

    # now take first derivate of event rates
    # first take t>0
    t_not_index = -1
    for i, time in enumerate(time_bins_x_axis):
        if time > 0:
            t_not_index = i
            break

    rates = all_plot_data[t_not_index:]
    # search for maximum
    max_rate_i = 0
    max_rate = 0
    for i, r in enumerate(rates):
        if r > max_rate:
            max_rate_i = i
            max_rate = r
    # linearly-spaced bins here, so can just double the index to approximate the burst window


    # t.create_regular_plot(all_plot_data,
    #                       handlers.same_axes(),
    #                       f'*Detectors ER {config.model_type} {config.transformation} {config.model_file_paths[number].split("/")[-1]}.png',
    #                       x_axis=time_bins_x_axis,
    #                       ylab='Event rate',
    #                       show=False
    #                       )
    #
    # # now renormalize and convert all points back to tuples
    # normalized = []
    # for point in all_plot_data:
    #     a = point[0]
    #     b = point[1]
    #     c = point[2]
    #     tot = a + b + c
    #     normalized.append((100 * a / tot, 100 * b / tot, 100 * c / tot))
    # all_plot_data = [tuple(point[0]) for point in all_plot_data]
    # t.create_regular_plot(normalized, handlers.same_axes(), f'{config.model_type} Super Normalized Ternary Points', 'Event Rate',
    #                       show=show_charts)

# we want to iterate over all snewpy models and their sets
# first find nmo data points
NMO: [float] = []
IMO: [float] = []
aggregate_detector(t.MetaAnalysisConfig(snewpy_models['Nakazato_2013'],[0],'NoTransformation'),0)


# for model in snewpy_models.keys():
#     for set in snewpy_models[model].file_paths:
#         # compute all event rates, then select based off of burst window
#