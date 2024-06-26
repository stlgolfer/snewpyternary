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
import pandas as pd

import snowglobes_wrapper
from matplotlib.widgets import Slider, Button
import matplotlib.animation as animation

sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import shared_plotting_script,generate_heatmap_dict,consolidate_heatmap_data
from model_wrappers import snewpy_models, sn_model_default_time_step
import click
from IPython import display

#simulation details
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'

# pulled the following from snowglobes-snewpy integration code
flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
complete_transform_list = list(flavor_transformation_dict.keys())
transforms_to_analyze = complete_transform_list # ['NoTransformation'] #['AdiabaticMSW_NMO','AdiabaticMSW_IMO','NoTransformation']

# global params
show_charts: bool = True
use_log: bool = True
use_cache: bool = True
use_presn: bool = False
use_all_submodules: bool = False
d: int = 10 # in pc, distance to SN
use_heatmap: bool = False

_colors = ['RED', 'GREEN', 'BLUE']
_run_cave_parameters: [(str, float)] = []

# Nt constants
NT_NUX = 1.004E35
NT_NUE = 6.04E34
NT_ANUE = 6.691E35

def _column_sum_proxy(data: dict, channels: ([str], [str], [str])) -> [float]:
    '''
    Like _sum_proxy, but does a column sum instead

    Parameters
    ----------
    data
    channels

    Returns
    -------
    A 3xN (N corresponding to number of time bins) array of the energy-dependent counts

    '''
    proxies = [
        np.zeros_like(data['Energy']),
        np.zeros_like(data['Energy']),
        np.zeros_like(data['Energy'])
    ]
    for index, proxy_flavor in enumerate(list(channels)):
        # print(proxy_flavor)
        # proxy_flavor has type [str]
        if len(proxy_flavor) > 0:
            sum = np.zeros_like(data['Energy'])

            for c in proxy_flavor:
                sum = np.add(sum,data[c])
            proxies[index] = sum

    return proxies

def process_detector(config: t.MetaAnalysisConfig, set_no: int, detector: str) -> None:
    global use_log

    plot_data, raw_data, l_data = t.create_detector_event_scatter(
        config.model_file_paths[set_no],
        config.model_type,
        detector,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transformation=config.transformation,
        data_calc=config.proxyconfig.build_detector_profiles()[detector]['chans_to_add'],
        use_cache=use_cache,
        log_bins=use_log,
        presn=use_presn,
        smearing='smeared'
    )

    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[set_no],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=use_log,
        presn=use_presn
    )

    # 3d spectra in time plot (https://matplotlib.org/stable/gallery/mplot3d/bars3d.html#sphx-glr-gallery-mplot3d-bars3d-py)
    # going to calculate here now

    spt_fig = plt.figure()
    spt_ax = spt_fig.add_subplot()

    sptc_index = {
        'ar40kt': {
            'index': 1,
            'proxy': 'nue'
        },
        'wc100kt30prct': {
            'index': 2,
            'proxy': 'anue'
        },
        'scint20kt': {
            'index': 0, # 0 for nux as usual
            'proxy': 'nux'
        }
    }

    # spt_title = f'{config.model_type} {detector}\n {str(config.proxyconfig)} {config.transformation} {"Logged" if use_log else "Linear"} Bins {sptc_index[detector]["proxy"]} Spectra{" PreSN" if use_presn else ""}'
    spt_title = f'{config.stringify(submodel=set_no)} {detector}'
    spt_full_content = []
    # now go through the l_data, which has rows containing dict_data
    with tqdm(total = len(l_data)) as pbar:
        for time_bin_no, dict_data in enumerate(l_data):
            spt_content = _column_sum_proxy(
                    dict_data,
                    config.proxyconfig.build_detector_profiles()[detector]['chans_to_add']
            )

            if (time_bin_no == 0):
                spt_full_content = spt_content[sptc_index[detector]['index']]
            else:
                spt_full_content = np.column_stack((spt_full_content, spt_content[sptc_index[detector]['index']]))
            pbar.update(1)
    spt_ax.set_ylabel('Energy (GeV)')
    spt_ax.set_xlabel(r'Time + $t_0$ (s)')
    # spt_ax.set_zlabel('Event rate')
    spt_ax.set_title(spt_title)
    # x dim should be energy bins, y should be time?
    __X, __Y = np.meshgrid((time_bins_x_axis/u.s), l_data[0]['Energy'])

    pc = spt_ax.contourf(__X, __Y, spt_full_content,10)
    spt_fig.colorbar(pc, shrink=0.75, location='right', label='Event Count/(GeV*s)', format='%.0e')
    spt_ax.set_xlim(0.0001, 20)
    spt_ax.set_xscale('log')

    # dump the figure for later arrangement
    spt_fig.savefig(f'./spectra/{t.clean_newline(spt_title)}.png')
    # pickle.dump(spt_fig, open(f'./spectra/{t.clean_newline(spt_title)}.pickle', 'wb'))

    if show_charts:
        spt_fig.show()

    figure, tax = t.create_default_detector_plot(
        raw_data,
        config.proxyconfig.build_detector_profiles()[detector]['axes'](),
        f'{config.model_type} {detector}\n{str(config.proxyconfig)} {config.transformation} {"Logged" if use_log else "Linear"} Bins Ternary{" PreSN" if use_presn else ""}',
        show=show_charts,
        save=True
    )

    # yeah, we're going to reprocess the flux because inefficient code is the best code for the unfolding
    # N_det is contained in the raw_data points
    # so, we should just be able to do a nice numpy element-wise operation?
    # but remember that the flux (1,2,3) has a different ordering, so we'll need to reshuffle
    # so just do it for water for now

    detector_to_index = {
        'wc100kt30prct': {
            'Ndet': 2,
            'phi_t': 1,
            'Nt': config.proxyconfig.Nt_wc100kt30prct()[2],#2,
            'proxy_name': 'IBD'
        },
        'scint20kt':{
            'Ndet': 0,
            'phi_t':0,
            'Nt': config.proxyconfig.Nt_scint20kt()[0],# 0,
            'proxy_name': 'NC'
        },
        'ar40kt': {
            'Ndet': 1,
            'phi_t':2,
            'Nt': config.proxyconfig.Nt_ar40kt()[1],#1,
            'proxy_name': 'Ar40'
        }
    }
    flux_scatter, flux_raw, flux_l_data = process_flux(config, set_no)
    N_det = np.transpose(np.array(raw_data))[detector_to_index[detector]['Ndet']] # this is the summed values for each time slice
    phi_t = np.transpose(np.array(flux_raw))[detector_to_index[detector]['phi_t']]
    n_targets = detector_to_index[detector]['Nt']
    # ok now I see the problem: the flux and the actual detector data are binned differently
    # the flux data has smaller energy bin width so the dE_v isn't the same
    flux_energy_spectra = np.linspace(0, 100, 501) #* MeV  # 1MeV

    #region anue spectrogram
    if detector == 'wc100kt30prct':
        # make a spectrogram of the flux for just anue
        flux_spect_fig, flux_spect_ax = plt.subplots(1, 1)
        flux_spectrogram = flux_l_data[0][4]
        for flux_spect_anue_bin in flux_l_data[1:]:
            flux_spectrogram = np.column_stack((flux_spectrogram, flux_spect_anue_bin[4]))

        flux_spect_ax.set_ylabel('Energy (MeV)')
        flux_spect_ax.set_xlabel('Time (s)')
        flux_spect_ax.set_title(r'$\bar{\nu_e}$ Flux Spectrogram')
        # x dim should be energy bins, y should be time?
        __X, __Y = np.meshgrid((time_bins_x_axis / u.s), flux_energy_spectra)

        flux_spect_ax.set_xlim(0.0001, 20)
        # flux_spect_ax.set_ylim(5, 100)
        # pcolormesh
        flux_spect_pc = flux_spect_ax.contourf(__X, __Y, flux_spectrogram,50) # , cmap=plt.cm.get_cmap('binary')
        # flux_spect_pc.set_clim(vmin=0,vmax=10000000.0)
        # going to add a slider for the max value
        # z_axes_max_slider = Slider(slider_axes, 'Blue', 0, 1e9)
        def __update(frame_num):
            flux_spect_pc.set_clim(vmin=0,vmax=frame_num)
        # z_axes_max_slider.on_changed(__update)
        flux_spect_fig.colorbar(flux_spect_pc, shrink=0.75, location='right', label=r'Neutrinos/(${cm}^2$*MeV*s)', format='%.0e')
        flux_spect_ax.set_xscale('log')

        # going to try an animation
        ani = animation.FuncAnimation(flux_spect_fig, __update,np.linspace(2e5,7e7,200), repeat=True)
        # code courtesy of Josh Q.
        print('Saving spectrogram video animation...')
        # writermp4 = animation.FFMpegWriter(fps=5)
        # ani.save('video.mp4', writer=writermp4)
        print('...Done')

        flux_spect_fig.savefig(f'./spectra/{config.stringify(submodel=set_no)} anue flux spectrogram.png')
        # pickle.dump(flux_spect_fig, open(f'./spectra/{config.stringify(submodel=set_no)} anue flux spectrogram.pickle', 'wb'))

        if show_charts:
            flux_spect_fig.show()
    #endregion


    return plot_data, raw_data, l_data

def process_flux(config: t.MetaAnalysisConfig, set_no: int, divide=1):

    flux_scatter_data, raw_data, labeled = t.create_flux_scatter(
        config.model_file_paths[set_no],
        config.model_type,
        config.model,
        deltat=sn_model_default_time_step(config.model_type),
        transform=config.transformation,
        use_cache=use_cache,
        log_bins=use_log,
        presn=use_presn
    )

    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[set_no],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=use_log,
        presn=use_presn
    )

    plt.figure()
    time_axis = np.arange(1, len(time_bins_x_axis) + 1, step=1)
    plt.plot(time_axis, time_bins_x_axis)
    plt.xlabel('Time Bin No')
    plt.ylabel('Calculated Time Coordinate')
    time_bin_plot_title = f'{config.model_type} {config.model_file_paths[set_no].split("/")[-1]} Truth Flux {"Logged" if use_log else "Linear"} Time Bins {config.transformation}{" PreSN" if use_presn else ""}.png'
    plt.title(time_bin_plot_title)
    plt.savefig(f'./plots/{time_bin_plot_title}')

    # now save this data
    axes_labels = data_handlers.ConfigBestChannel().flux_axes()
    data_reshuffled = np.vstack(raw_data).T
    df = pd.DataFrame({
        'time':time_bins_x_axis,
        f'raw_data_nux':data_reshuffled[0],
        f'raw_data_anue': data_reshuffled[1],
        f'raw_data_nue': data_reshuffled[2]
    })
    df.to_csv(f'./flux_saves/{config.stringify(set_no)}_flux_save.csv')
    return flux_scatter_data, raw_data, labeled

def remap_dict(dictionary,newval):
    # remaps a dictionary's 1 value to a different value
    new_dict = {}
    for k in dictionary.keys():
        if dictionary[k] == 1:
            new_dict[k] = newval
        else:
            new_dict[k] = 0
    return new_dict

def t_normalize(raw_data):
    normalized = np.zeros_like(raw_data)
    for point in range(len(normalized)):
        a = raw_data[point][0]
        b = raw_data[point][1]
        c = raw_data[point][2]
        tot = a + b + c
        normalized[point] = (100 * a / tot, 100 * b / tot, 100 * c / tot)
        #normalized.append()
    return normalized

def ternary_distance(p1: tuple, p2: tuple):
    '''
    Calculate the distance between two TERNARY distances
    Parameters
    ----------
    p1 a normalized t-space point
    p2 a normalized t-space point

    Returns
    -------

    '''
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)

def ternary_subtract(p1: tuple, p2: tuple):
    #p1 - p2
    return (p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2])

def ternary_dotproduct(p1: tuple, p2: tuple):
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2]

def aggregate_detector(config: t.MetaAnalysisConfig, number: int, colorid: int, tax: TernaryAxesSubplot, cum_sum_tax: TernaryAxesSubplot) -> None:
    # flux_scatter, flux_raw, flux_l_data = process_flux(config, number)

    # print out information of the set
    print(config.model(config.model_file_paths[number]))

    p_data, r_data, l_data = process_detector(config, number, 'ar40kt')
    # need to convert data to an array
    all_plot_data = [list(key) for key in r_data]  # going to take each detector and add them up

    for detector in ['wc100kt30prct', 'scint20kt']:
        p_data, r_data, l_data = process_detector(config, number, detector)
        all_plot_data = all_plot_data + np.asarray([list(key) for key in r_data])

    # now get the time bins
    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[number],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=use_log,
        presn=use_presn
    )
    # r'$\nu_x$ Proxy', r'$\nu_e$ Proxy', r'$\bar{\nu_e}$ Proxy'

    # region also create a cumulative plot
    nux_proxy_cumsum = np.cumsum(list(list(zip(*all_plot_data))[0]))
    nue_proxy_cumsum = np.cumsum(list(list(zip(*all_plot_data))[1]))
    anue_proxy_cumsum = np.cumsum(list(list(zip(*all_plot_data))[2]))

    # create a new ternary diagram for the cumsum
    cumsum_normalized = []
    for ci in range(len(nux_proxy_cumsum)):
        ci_total = nux_proxy_cumsum[ci] + nue_proxy_cumsum[ci] + anue_proxy_cumsum[ci]
        cumsum_normalized.append(
            (100*nux_proxy_cumsum[ci]/ci_total, 100*nue_proxy_cumsum[ci]/ci_total, 100*anue_proxy_cumsum[ci]/ci_total))
    # endregion

    # now renormalize and convert all points back to tuples
    normalized = t_normalize(all_plot_data)

    # going to try dynamically sized points between lines?
    widths = np.linspace(0.01, 1, num=len(normalized))
    cs_widths = np.linspace(0.01,1,num=len(cumsum_normalized))
    for p in range(len(normalized) - 1):
        if (p + 1 >= len(normalized)):
            break
        tax.line(normalized[p], normalized[p + 1], color=(widths[p] if colorid == 0 else 0, widths[p] if colorid == 1 else 0, widths[p] if colorid == 2 else 0, 1), linestyle=':', linewidth=3)
        # TODO: fix cumsum
        # cum_sum_tax.line(cumsum_normalized[p], cumsum_normalized[p + 1], color=(
        # cs_widths[p] if colorid == 0 else 0, cs_widths[p] if colorid == 1 else 0, cs_widths[p] if colorid == 2 else 0, 1),
        #          linestyle=':', linewidth=3)
    # tax.scatter(normalized, color='blue')

    if use_heatmap:
        print('Calculating errorbar heatmap...')
        # more information on colormaps can be found here:
        # https://matplotlib.org/stable/tutorials/colors/colormaps.html#diverging
        tax.heatmap(generate_heatmap_dict(all_plot_data,t_normalize(all_plot_data)), cmap=plt.get_cmap('PiYG'))
        print('...Done')

    # here we also want to calculate the cave parameter
    # first need the cs track length, which is a ternary-space line integral. yikes
    line_integral = 0
    for bin_no, p in enumerate(cumsum_normalized[0:-1]):
        # skip last element since we're doing a delta
        line_integral = line_integral + ternary_distance(p, cumsum_normalized[bin_no+1])
    # now measure the cave parameter
    cave_param = ternary_distance(cumsum_normalized[0], cumsum_normalized[-1])/line_integral
    global _run_cave_parameters
    _run_cave_parameters.append((config.model_file_paths[number].split("/")[-1], cave_param))

def process_transformation(config: t.MetaAnalysisConfig):
    print(f'Now processing {config.model_type}')
    scale=100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale/10)
    title=t.clean_newline(f'{config.model_type} *Detectors Folded {config.transformation} {str(config.proxyconfig)}\n {"Logged" if use_log else "Linear"} Bins{" PreSN" if use_presn else ""}{" AS" if use_all_submodules else ""} Ternary')
    tax.set_title(title)
    # data is organized in top, right, left
    # apparently this is in flux formatting
    tax.bottom_axis_label('nux')
    tax.right_axis_label('nuebar')
    tax.left_axis_label('nue')

    # now draw the cumsum chart
    cumsum_figure, cum_sum_tax = ternary.figure(scale=100)
    cum_sum_tax.boundary(linewidth=2.0)
    cum_sum_tax.gridlines(color="blue", multiple=100 / 10)
    cumsum_title = f'{config.model_type} *Detectors Folded Cumsum {config.transformation} {str(config.proxyconfig)}\n {"Logged" if use_log else "Linear"} Bins{" PreSN" if use_presn else ""}{" AS" if use_all_submodules else ""} Ternary'
    cum_sum_tax.set_title(cumsum_title)
    # data is organized in top, right, left

    cum_sum_tax.bottom_axis_label('nux')
    cum_sum_tax.right_axis_label('nue')
    cum_sum_tax.left_axis_label('nuebar')

    # heatmap line goes here
    #timemap = {}
    #colorspace = [c for c in range(len(all_plot_data))]
    #for j, n in enumerate(normalized):
        #timemap[n] = float(j+10)
        # then vary each coordinate
        #p_i = n # initial tuple coordinate in ternary space
        #for r in range(heatmap_radius):
            #p_pos = (int(p_i[0]) + r,int(p_i[1]) + r, int(p_i[2]) +r)
            #p_neg = (int(p_i[0]) - r,int(p_i[1]) - r, int(p_i[2]) -r)
            # add to timemap
            #timemap[p_pos] = j
            #timemap[p_neg] = j
            
    #if show_time_heatmap == True:
    #    tax.heatmap(timemap)

    # need to create different colors if not using all submodels
    if use_all_submodules:
        for set_number in config.set_numbers:
            aggregate_detector(config, set_number, 0, tax, cum_sum_tax)
    else:
        colorid: int = 0
        f = open(f'./all_detector_plots/{title}.txt', 'w')
        for set_number in config.set_numbers:
            # write model metadata as well
            f.write(f'{_colors[colorid]}\n==========\n')
            f.write(repr(config.model(config.model_file_paths[set_number])))
            f.write('\n\n')
            aggregate_detector(config, set_number, colorid, tax, cum_sum_tax)
            colorid += 1
        f.close()

    #tax.scatter(points=normalized)
    
    tax.ticks(axis='lbr', linewidth=1, multiple=scale/10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off') # disables regular matlab plot axes
    tax.savefig(f'./all_detector_plots/{title}')
    
    if show_charts == True:
        figure.show()

    # region show cum sum ternary
    cum_sum_tax.ticks(axis='lbr', linewidth=1, multiple=100 / 10)
    cum_sum_tax.clear_matplotlib_ticks()
    cum_sum_tax.get_axes().axis('off')  # disables regular matlab plot axes
    # cum_sum_tax.savefig(f'./all_detector_plots/{t.clean_newline(cumsum_title)}')

    # if show_charts == True:
        # cumsum_figure.show()
    #endregion

    #region cave parameter output file
    # if do_unfold:
    #     global _run_cave_parameters
    #     # now write cave params to file only if we are doing unfolding. otherwise it doesn't make much sense
    #     cp_title = t.clean_newline(
    #         f'{config.model_type} *Detectors Cave Parameter {config.transformation} {str(config.proxyconfig)}\n {"Logged" if use_log else "Linear"} Bins{" PreSN" if use_presn else ""}')
    #     cp_f = open(f'./cave_parameters/{cp_title}.txt', 'w')
    #     for cp_run in _run_cave_parameters:
    #         cp_f.write(f'{cp_run[0]},{cp_run[1]}\n')
    #     cp_f.close()
    #     # reset run cave parameters for next transformation
    #     _run_cave_parameters = []
    #endregion

# process_transformation(t.MetaAnalysisConfig(snewpy_models['Bollig_2016'], 'NoTransformation'))
@click.command()
@click.argument('models',required=True,type=str,nargs=-1)
@click.option('--showc',default=False,type=bool,help='Whether to show generated plots or not. Will always save and cache')
@click.option('-p',required=False, multiple=True, type=str, default=['NoTransformation'], help='Prescriptions to use')
@click.option('--distance',default=10,type=int,help='The distance (in kPc) to the progenitor source')
@click.option('--uselog',default=True,type=bool, help='Use logarithmically spaced (and shifted) time bins')
@click.option('--setno', required=False, default=[0],type=int,multiple=True, help='Model set index. See model_wrappers.py')
@click.option('--allsubmodels', required=False, type=bool, help="Cannot specify setno and this at the same time. Runs setnos available submodules")
@click.option('--cache', required=False, default=True, type=bool, help='If true, use cache')
@click.option('--presn', required=False, default=False, type=bool, help='If true, compute time bins from t<=0')
@click.option('--tflux', required=False, default=False, type=bool, help='If true, only calculate the truth flux. set numbers are not superimposed')
@click.option('--detproxy', required=False, type=str, default='AgDet', help='Detector proxy configuration. Options: AgDet or BstChnl')
@click.option('--heatmap', required=False, type=bool, default=False, help='Include heatmaps. Can only process one submodel at a time')
def start(showc,models,distance,uselog,p, setno, allsubmodels, cache, presn, tflux, detproxy, heatmap):
    if detproxy != 'BstChnl':
        raise ValueError('Use of detproxy option is currently deprecated. Only using BstChnl')

    global show_charts
    show_charts = showc
    
    global d
    d = distance
    
    global use_log
    use_log = uselog

    global use_cache
    use_cache = cache

    global use_presn
    use_presn = presn

    global use_all_submodules
    use_all_submodules = allsubmodels

    global use_heatmap
    use_heatmap = heatmap

    if use_heatmap and len(setno) > 1:
        raise ValueError('Can only process heatmap for one submodel at a time')

    # check set numbers
    if len(setno) > 3:
        raise ValueError("Can only superimpose a maximum of 3 sets onto one chart")

    if use_presn and use_log:
        raise RuntimeError("Using log time bins on presn calculations are not currently supported")

    # want to iterate by model, then prescription, then set

    for model in models:
        for prescription in (flavor_transformation_dict.keys() if model == "ALL" else p):
            # check to see if valid model set number

            for no in setno:
                if no >= len(snewpy_models[model].file_paths):
                    raise ValueError(f"Invalid model set id. Max is {len(snewpy_models[model].file_paths)-1}")

            # remove multithreading for now. run sequentially
            proxy = data_handlers.ConfigAggregateDetectors() if detproxy == 'AgDet' else data_handlers.ConfigBestChannel()
            # setno is a collection of indexes in the model_wrapper lookups. if all models are selected then linspace
            selected_setnos = [x for x in range(len(snewpy_models[model].file_paths))] if allsubmodels else setno
            config = t.MetaAnalysisConfig(snewpy_models[model], selected_setnos, prescription,proxy_config=proxy)

            if tflux:
                for num in selected_setnos:
                    process_flux(config, num)
            else:
                process_transformation(config)

            # proc = mp.Process(target=process_transformation, args=[])
            # proc.start()
            # proc.join()
            
if __name__ == '__main__': # for multiprocessing
    start()

# for d in handlers.supported_detectors:
# for d in handlers.supported_detectors:
#     nmo_figure, nmo_tax, nmo_plot_data, nmo_raw_data, nmo_hp = process_detector(d,transforms_to_analyze[0])
#     imo_figure, imo_tax, imo_plot_data, imo_raw_data, imo_hp = process_detector(d,transforms_to_analyze[1])
    
#     # just gonna have to create a new graph for combining things
#     figure, tax = ternary.figure(scale=100)
#     tax.boundary(linewidth=2.0)
#     tax.gridlines(color="black", multiple=10)
#     title = f'{model_type} {d} NMO vs IMO'
#     tax.set_title(title)
#     # data is organized in top, right, left
    
#     axes_titles = profiles[d]['axes']()
#     tax.bottom_axis_label(axes_titles[0])
#     tax.right_axis_label(axes_titles[1])
#     tax.left_axis_label(axes_titles[2])
#     d1 = nmo_hp.copy()
#     d2 = imo_hp.copy()
#     newval = 0.6
#     d2 = remap_dict(d2,newval)
#     tax.heatmap(consolidate_heatmap_data(d1,d2),cbarlabel='In 68% CL')
    
#     tax.scatter(points=nmo_plot_data, color='red')
#     tax.scatter(points=imo_plot_data, color='blue')
    
#     tax.ticks(axis='lbr', linewidth=1, multiple=10)
#     tax.clear_matplotlib_ticks()
#     # tax.get_axes().axis('off') # disables regular matlab plot axes
#     tax.show()
#     tax.savefig(f'./plots/{title}.png')

# for trans in transforms_to_analyze: process_flux(trans)

# if multithreading==True:
#     for detector in profiles.keys():
#         proc = Process(target=process_detector, args=[detector])
#         proc.start()
#         proc.join()
# else:
#     for detector in detectors.keys():
#         process_detector(detector)
