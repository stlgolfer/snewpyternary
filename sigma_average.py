import click

import snewpyternary as t
import snowglobes_wrapper
from model_wrappers import snewpy_models, sn_model_default_time_step
import data_handlers
from meta_analysis import process_flux, process_detector, t_normalize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
from astropy import units as u
from tqdm import tqdm

def estimate_cxn(
        config: t.MetaAnalysisConfig,
        flux_chan: int,
        cxn_truth_fname: str,
        cxn_truth_chan_key: str,
        det_name: str,
        det_chan_name: str,
        Nt: float
):
    submodel_number = config.set_numbers[0]
    time_bins_x_axis, dt_not_needed = snowglobes_wrapper.calculate_time_bins(
        config.model_file_paths[submodel_number],
        config.model_type,
        deltat=sn_model_default_time_step(config.model_type),
        log_bins=True,
        presn=False
    )

    times_unitless = [m.value for m in time_bins_x_axis]
    temp1 = np.array([times_unitless[0] + 0.5])
    temp2 = times_unitless[1:]
    temp3 = times_unitless[0:-1]
    temp4 = np.subtract(np.array(temp2), np.array(temp3))
    dts = np.concatenate((np.array(temp1), temp4))

    flux_scatter_data, raw_data, labeled = process_flux(config, submodel_number)

    # in a time bin, the truth flux is binned in GeV
    labeled_transposed = np.transpose(labeled)
    plot_data, raw_data_det, l_data = process_detector(config, submodel_number, det_name)

    #region make spectrograms and save them in the right spot
    spt_full_content = []
    # now go through the l_data, which has rows containing dict_data

    #endregion

    #region flux ternary diagram
    flux_td_title = f'{config.stringify(config.set_numbers[0])} Truth Flux'
    # need to divide the nux data by 6
    flux_td_fig, flux_td_tax = t.create_default_flux_plot(t_normalize([(x[0]/6,x[1],x[2]) for x in raw_data]), flux_td_title)
    flux_td_tax.show()
    flux_td_tax.savefig(f'./plots/{flux_td_title} Ternary Diagram.png')
    #endregion

    # region let's just try to plot phi_t over time, flux vs time, and flux vs energy
    phi_t_fig, phi_t_ax = plt.subplots(1, 1)

    phi_t_over_time = np.zeros_like(times_unitless)
    flux_vs_time = np.zeros_like(times_unitless)
    for l in range(len(labeled)):
        nux_fluence = labeled[l][1] + labeled[l][2] + labeled[l][3] + labeled[l][4] + labeled[l][5] + \
                      labeled[l][6]
        phi_t_over_time[l] = np.sum(labeled[l][flux_chan] if det_name != "scint20kt" else nux_fluence) * 0.2e-3  * dts[l]
        flux_vs_time[l] = np.sum(labeled[l][flux_chan] if det_name != "scint20kt" else nux_fluence) * 0.2e-3

    if det_name == 'scint20kt':
        warnings.warn("scint20kt uses a slightly different calculation, so don't be surprised if the AUC is wrong")

    phi_t_ax.scatter(times_unitless, phi_t_over_time)
    phi_t_ax.set_title(rf'$\phi_t$ for {config.model_type} {config.set_numbers[0]} Fluence')
    phi_t_ax.set_xlabel(r"Time + $t_0$ (s)")
    phi_t_ax.set_ylabel(r'$neutrinos/cm^2$')
    phi_t_ax.set_xscale('log')
    phi_t_inset_ax = phi_t_ax.inset_axes([0.2,0.2,0.5,0.5])
    phi_t_inset_ax.scatter(times_unitless, phi_t_over_time)
    phi_t_inset_ax.set_xscale('log')
    phi_t_inset_ax.set_yscale('log')
    phi_t_fig.tight_layout()
    phi_t_fig.show()

    flux_vs_time_fig, flux_vs_time_ax = plt.subplots(1, 1)
    flux_vs_time_ax.bar(times_unitless, flux_vs_time, width=dts)
    flux_vs_time_ax.set_title("Flux vs time, integrated over energy")
    flux_vs_time_ax.set_xlabel("Time Mid-Point (s)")
    flux_vs_time_ax.set_ylabel(r'$\frac{neutrinos}{{cm^2}*s}$')
    flux_vs_time_fig.tight_layout()
    flux_vs_time_AUC = np.sum(np.multiply(flux_vs_time, dts))
    print(f"AUC for Flux vs time is {flux_vs_time_AUC}")
    flux_vs_time_fig.show()

    flux_vs_energy = np.zeros_like(labeled[0][0])
    for en in range(len(flux_vs_energy)):
        flux_vs_energy[en] = np.sum(np.multiply(labeled_transposed[en][flux_chan], dts))
    flux_vs_energy_fig, flux_vs_energy_ax = plt.subplots(1, 1)
    flux_vs_energy_ax.bar(labeled[0][0], flux_vs_energy, width=0.2e-3)
    flux_vs_energy_ax.set_title("Flux vs energy, integrated over time")
    flux_vs_energy_ax.set_xlabel('Mid-Point Energy (GeV)')
    flux_vs_energy_ax.set_ylabel(r'$\frac{neutrinos}{{cm}^2 * GeV}$')
    # then find the AUC
    flux_vs_energy_AUC = np.sum(flux_vs_energy) * 0.2e-3
    print(f"AUC for Flux vs energy is {flux_vs_energy_AUC}")
    print(f'AUC for phi_t is {np.sum(phi_t_over_time)}')
    flux_vs_energy_fig.show()

    print(f"AUC ratio fve/fvt = {flux_vs_energy_AUC / flux_vs_time_AUC}")
    # endregion

    # truth_flux_fig, truth_flux_ax = plt.subplots(1,1)
    # truth_flux_ax.scatter()

    # make a spectrogram of the actual xscn from snowglobes (it will be time-invariant)
    cxn_truth = pd.read_csv(f'./snowglobes_cxns/{cxn_truth_fname}.csv')  # only want the 'nu_e_bar' channel
    # now we'll also need to downsample the cross since the cross section has far more bins (around double) than flux
    # this is also time-invariant so only need to compute once
    print("Loaded truth xscn")

    cxn_truth_energy = 10 ** np.array(cxn_truth['energy'])

    truth_calculation = np.multiply(
        cxn_truth_energy,
        np.array(cxn_truth[cxn_truth_chan_key])
    )

    truth_calculation_no_downsample = np.copy(truth_calculation) * 1e-38

    truth_calculation = np.nan_to_num(truth_calculation, nan=0.0)

    # region reconstruct cxn from Ndet to ensure correct smearing matrix was used
    # pick 0th bin
    nux_fluence = labeled[0][1] + labeled[0][2] + labeled[0][3] + labeled[0][4] + labeled[0][5] + labeled[0][6]
    flux_interpolated = np.interp(
        l_data[0]['Energy'],
        labeled[0][0],
        nux_fluence if det_name == 'scint20kt' else labeled[0][flux_chan]
    )

    cxn_reconstructed = np.divide(
        l_data[0][det_chan_name],
        flux_interpolated
    ) / (Nt * 2.5)
    recon_cxn_fig, recon_cxn_ax = plt.subplots(1,1)
    recon_cxn_ax.scatter(l_data[0]['Energy'], cxn_reconstructed)
    recon_cxn_ax.set_xscale('log')
    recon_cxn_ax.set_yscale('log')
    recon_cxn_fig.show()
    # endregion

    # downsample the truth flux with 1e38 cm^2 scale since python might have a hard time interpolating such small values
    truth_calculation = np.interp(labeled[0][0], cxn_truth_energy, truth_calculation) * 1e-38
    print("Truth CXN Downsampled")

    # region plot nu_e_bar_cxn_truth
    cxn_truth_fig, cxn_truth_ax = plt.subplots(1, 1)
    cxn_truth_ax.scatter(cxn_truth_energy, truth_calculation_no_downsample, label="CXN No Downsample")
    cxn_truth_ax.scatter(labeled[0][0], truth_calculation, label="CXN Downsampled", alpha=0.4)
    cxn_truth_ax.set_xlabel('Energy (GeV)')
    cxn_truth_ax.set_ylabel(r'${cm}^2$')
    cxn_truth_ax.set_title(r'$\sigma$ Truth')
    cxn_truth_ax.set_xlim(1e-3, 1e-1)
    cxn_truth_ax.set_ylim(1e-44, 1e-39)
    cxn_truth_ax.set_xscale('log')
    cxn_truth_ax.set_yscale('log')
    cxn_truth_fig.legend()
    cxn_truth_fig.show()
    # endregion

    # region plot truth flux average energy over time
    truth_flux_average_energy = np.zeros_like(times_unitless)
    for lt in range(len(labeled)):
        # find int E*F(E) dE / int F(E) dE, though in all fairness dE will cancel here
        truth_flux_average_energy[lt] = np.sum(np.multiply(labeled[lt][0], labeled[lt][flux_chan] * dts[lt])) / np.sum(
            labeled[lt][flux_chan] * dts[lt])

    tf_average_energy_fig, tf_average_energy_ax = plt.subplots(1, 1, figsize=(8, 8))
    tf_average_energy_ax.scatter(times_unitless, truth_flux_average_energy)
    tf_average_energy_ax.set_xscale('log')
    tf_average_energy_ax.set_xlabel('Time (s)')
    tf_average_energy_ax.set_ylabel('Flux-Weighted Average Energy (GeV)')
    tf_average_energy_ax.set_title(r'$\mathbb{E}(E_\nu|F(E_\nu))$')
    tf_average_energy_fig.show()
    # endregion

    # we want to make the flux average cross-section by using the cross-section from snowglobes
    # go through each time bin

    sigma_average = np.zeros_like(times_unitless)
    for time in range(len(times_unitless)):
        nux_fluence = labeled[time][1] + labeled[time][2] + labeled[time][3] + labeled[time][4] + labeled[time][5] + \
                      labeled[time][6]
        sigma_average[time] = np.sum(
            np.multiply(
                dts[time] * 0.2e-3 * (nux_fluence if det_name == 'scint20kt' else labeled[time][flux_chan]),
                truth_calculation)
        ) / phi_t_over_time[time]

    # now plot the flux-averaged cxn
    fig, ax = plt.subplots(1, 1)
    ax.scatter(times_unitless, sigma_average)
    ax.set_ylabel(r'$<\sigma>$')
    ax.set_xlabel('Time (s)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(rf'$<\sigma>$ for Nakazato {config.set_numbers[0]} {cxn_truth_chan_key}')
    fig.show()

    #region unfold attempt
    # at this point, we have everything we need to also unfold. the form should be similar to the phi_t calculation
    # for now, just try unfolding anue
    phi_est_unfolded = np.zeros_like(times_unitless)
    Ndet_over_time = np.zeros_like(times_unitless)
    for phi_est_time_bin in range(len(times_unitless)):
        if det_name == 'scint20kt':
            l_data[phi_est_time_bin][det_chan_name][l_data[0]['Energy']*1000 < 15] = 0

        Ndet = np.sum(
            l_data[phi_est_time_bin][det_chan_name]
        ) # dts[phi_est_time_bin] *
        #TODO: need to figure out if this needs to be commented out or not. it looks like rishi's code
        # considers when they're maybe still binned? no the integration part should be right
        # look at line 220 the integration already happens so this should in fact be total count and units line up
        Ndet_over_time[phi_est_time_bin] = Ndet #/(0.2e-3 * dts[phi_est_time_bin])

        phi_est_unfolded[phi_est_time_bin] = Ndet / (Nt * sigma_average[phi_est_time_bin])
        # ok double check where the 0.2e-3 * dts are happening because i bet something isn't getting canceled somewhere
    unfold_fig, unfold_ax = plt.subplots(1,1)
    unfold_ax.set_xscale('log')
    unfold_ax.scatter(times_unitless, phi_est_unfolded, label='Unfolded')
    unfold_ax.scatter(times_unitless, phi_t_over_time, label=r'$\phi_t$ Truth', marker='1', color='red', sizes=100*np.ones_like(phi_t_over_time)) #, alpha=0.2
    unfold_ax.set_title(f'{config.stringify(submodel_number)} {cxn_truth_chan_key} Unfolded')
    unfold_ax.legend()
    unfold_fig.show()
    #endregion

    #region also plot the Ndet_over_time as a general check
    ndet_fig, ndet_ax = plt.subplots(1)
    ndet_ax.set_xscale('log')
    ndet_ax.scatter(times_unitless, Ndet_over_time)
    ndet_ax.set_title(f'{config.stringify(submodel_number)} {cxn_truth_chan_key} Detector Count')
    ndet_ax.set_xlabel('Time + t0 (s)')
    ndet_ax.set_ylabel('Detector Count')
    ndet_fig.show()
    #endregion

    if __name__ == '__main__':
        # only save the sigmas in a csv if this is being run by itself
        print('Storing average sigma in ./sigmas...')
        df = pd.DataFrame()
        df['time'] = np.array(times_unitless) + 0.0001 + abs(times_unitless[0]) #TODO: shift axis to include more of the time window
        df['dt'] = dts
        df['sigma average'] = sigma_average
        df['unfolded'] = phi_est_unfolded
        df['Ndet'] = Ndet_over_time

        df.to_csv(f'./sigmas/{config.stringify(config.set_numbers[0])}_{cxn_truth_chan_key}_sigma_average.csv')
        print('Done')

    return sigma_average


@click.command()
@click.argument('model', required=True, type=str, nargs=1)
@click.option('-p', required=True, type=str)
@click.option('-flavor', required=True, type=str) # help='nue, anue, or nux'
@click.option('--submodel', required=False, type=int, help='Internal submodel index number', default=0)
def configure(model, p, flavor, submodel):
    warnings.warn('Will use global settings found in meta_analysis.py')
    warnings.warn('Only BstChnl configuration supported')
    '''
    Internal helper function to convert click parameters to config 
    Returns
    -------

    '''
    config = t.MetaAnalysisConfig(
        snewpy_models[model],
        [submodel],
        p,
        proxy_config=data_handlers.ConfigBestChannel()
    )

    flavor_to_config = {
        'nux': {
            'flux_chan': 1,
            'cxn_truth_fname': 'xs_nc_numu_C12',
            'cxn_truth_chan_key': 'nu_mu',
            'det_name': 'scint20kt',
            'det_chan_name': 'nc',
            'Nt': config.proxyconfig.Nt_scint20kt()[0]
        },
        'nue': {
            'flux_chan': 1,
            'cxn_truth_fname': 'xs_nue_Ar40',
            'cxn_truth_chan_key': 'nu_e',
            'det_name': 'ar40kt',
            'det_chan_name': 'nue_Ar40',
            'Nt': config.proxyconfig.Nt_ar40kt()[1]
        },
        'anue': {
            'flux_chan': 4,
            'cxn_truth_fname': 'xs_ibd',
            'cxn_truth_chan_key': 'nu_e_bar',
            'det_name': 'wc100kt30prct',
            'det_chan_name': 'ibd',
            'Nt': config.proxyconfig.Nt_wc100kt30prct()[2]
        }
    }

    estimate_cxn(
        config,
        flavor_to_config[flavor]['flux_chan'],
        flavor_to_config[flavor]['cxn_truth_fname'],
        flavor_to_config[flavor]['cxn_truth_chan_key'],
        flavor_to_config[flavor]['det_name'],
        flavor_to_config[flavor]['det_chan_name'],
        flavor_to_config[flavor]['Nt']
    )

if __name__ == '__main__':
    configure()
