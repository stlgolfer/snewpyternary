import snewpyternary as t
import snowglobes_wrapper
from model_wrappers import snewpy_models, sn_model_default_time_step
import data_handlers
from meta_analysis import process_flux, process_detector
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u

if __name__ == '__main__':
    config = t.MetaAnalysisConfig(
        snewpy_models['Nakazato_2013'],
        [0],
        'AdiabaticMSW_IMO',
        proxy_config=data_handlers.ConfigBestChannel()
    )

    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[0],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=True,
        presn=False
    )

    flux_scatter_data, raw_data, labeled = process_flux(config,0)
    # in a time bin, the truth flux is binned in GeV
    plot_data, raw_data_det, l_data, zeta, sigma_average = process_detector(config,0,'wc100kt30prct')

    #TODO: I suspect that when the multiplication between the cross-section and the truth flux happens
    # the energy units are not the same

    # truth_flux_fig, truth_flux_ax = plt.subplots(1,1)
    # truth_flux_ax.scatter()

    # make a spectrogram of the actual xscn from snowglobes (it will be time-invariant)
    ibd_cxn_truth = pd.read_csv('./time_cuts/xs_ibd.csv')  # only want the 'nu_e_bar' channel
    # now we'll also need to downsample the cross since the cross section has far more bins (around double) than flux
    # this is also time-invariant so only need to compute once
    print("Loaded truth xscn")

    nu_e_bar_cxn_truth = np.array(ibd_cxn_truth['nu_e_bar'])*10**np.array(ibd_cxn_truth['energy'])*1e-38
    ibd_cxn_truth_downsampled = np.interp(
        labeled[0][0],
        10**np.array(ibd_cxn_truth['energy']),
        nu_e_bar_cxn_truth
    )
    ibd_cxn_truth_downsampled = np.nan_to_num(ibd_cxn_truth_downsampled, nan=0.0)

    # make a plot to compare cxn truths to make sure it downsampled correctly
    cxn_truth_comp_fig, cxn_truth_downsampled_comp = plt.subplots(1, 1, figsize=(8, 8))
    cxn_truth_downsampled_comp.scatter(10 ** np.array(ibd_cxn_truth['energy']), nu_e_bar_cxn_truth) # this in GeV
    cxn_truth_downsampled_comp.scatter(labeled[0][0], ibd_cxn_truth_downsampled)
    cxn_truth_comp_fig.show()

    # we want to make the flux average cross-section by using the cross-section from snowglobes
    # go through each time bin

    # ok for now, let's just swap the phi_t with impulse 1s

    times_unitless = [m.value for m in time_bins_x_axis]
    temp1 = np.array([times_unitless[0] + 0.5])
    temp2 = times_unitless[1:]
    temp3 = times_unitless[0:-1]
    temp4 = np.subtract(np.array(temp2), np.array(temp3))
    dts = np.concatenate((np.array(temp1), temp4))

    sigma_average = np.zeros_like(times_unitless)
    conversion = ((0.095e-3)/0.0002)**-1
    for time in range(len(times_unitless)):
        # channel 4 in labeled is anue
        # phi_t = np.sum(labeled[time][4])*dts[time]*.0002 # 0.2 MeV energy bins for flux, but we're in GeV
        # what is the bin size of the sigma? -- it's 0.095 MeV
        impulse_substitute = np.zeros_like(labeled[time][4])
        impulse_substitute[100] = 1

        sigma_average[time] = np.sum(np.multiply(ibd_cxn_truth_downsampled, impulse_substitute))*(0.095) # / phi_t # should be times two for dEv

    # now plot the flux-averaged cxn
    fig, ax = plt.subplots(1,1)
    ax.scatter(times_unitless, sigma_average)
    ax.set_ylabel(r'$<\sigma>$')
    ax.set_xlabel('Time (s)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(rf'$<\sigma>$ for Nakazato 0 IBD for $\phi_t({labeled[0][0][100]} GeV) = 1$')
    fig.show()
