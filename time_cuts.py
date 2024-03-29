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
    plot_data, raw_data_det, l_data = process_detector(config,0,'scint20kt')

    # make a spectrogram of the actual xscn from snowglobes (it will be time-invariant)
    ibd_cxn_actual = pd.read_csv('./snowglobes_cxns/xs_nc_numu_C12.csv')  # only want the 'nu_e_bar' channel
    cxn_comp_fig, (cxn_truth_axes, cxn_actual_axes) = plt.subplots(1, 2, figsize=(16, 8))
    cxn_actual_spectrogram = None
    cxn_truth_spectrogram = None

    truth_calculation = np.multiply(
        10**np.array(ibd_cxn_actual['energy']),
        np.array(ibd_cxn_actual['nu_mu'])*1e-38
    )

    # want 0.05s, 0.07s, 0.1s, 0.3s, 0.5s, 1, 3, 5, 7, 10
    # time_cut_indexes = [114,117,121,135,142,153,170,178,184,190]
    for time in range(len(flux_scatter_data)):
        # fig, (flux_axes, ndet_axes, ratio_axes) = plt.subplots(1,3,figsize=(24,8))
        # want to superimpose flux and Ndet, so make sure GeV
        flux_energy_spectra = np.linspace(0, 100, 501)  # * MeV  # 1MeV
        # flux_axes.scatter(flux_energy_spectra,labeled[time][4],label='Flux')
        # flux_axes.set_title('Flux')
        # flux_axes.set_xlabel('Energy (MeV)')
        # flux_axes.set_ylabel('neutrinos/(cm^2 * MeV)')
        #
        # ndet_axes.scatter(l_data[time]['Energy']*1000,l_data[time]['ibd'], label='Ndet')
        # ndet_axes.set_title('Ndet')
        # ndet_axes.set_xlabel('Energy (MeV)')
        # ndet_axes.set_ylabel('Event count / MeV')

        # ndet_interpolated = np.interp(flux_energy_spectra,l_data[time]['Energy']*100,l_data[time]['ibd'])
        nux_fluence = labeled[time][1] + labeled[time][2] + labeled[time][3] + labeled[time][4] + labeled[time][5] + labeled[time][6]
        flux_interpolated = np.interp(l_data[time]['Energy']*1000, flux_energy_spectra, nux_fluence)
        cxn_reconstructed = np.divide(
            l_data[time]['nc'],
            flux_interpolated
        )/(config.proxyconfig.Nt_scint20kt()[0]*2.5) # divide by 2.5 since the energy bin widths are different

        # if we're using scint20kt, then we need to 0 out the cxn for energies less than 10 MeV (0.01 GeV)
        cxn_reconstructed[l_data[0]['Energy']*1000 < 15] = 0

        # ratio_axes.scatter(
        #     l_data[time]['Energy']*1000,
        #     cxn_reconstructed
        # )
        # ratio_axes.set_xscale('log')
        #
        # ratio_axes.set_title('Nt * Ndet/Flux')
        # ratio_axes.set_xlabel('Energy (MeV)')
        # ratio_axes.set_ylabel('Nt * Ndet/Flux')
        #
        # fig.suptitle(rf'Flux vs Ndet Spect. $t\approx{time_bins_x_axis[time]}$')
        #
        # fig.savefig(f'./time_cuts/{time_bins_x_axis[time]}.png')
        if time == 0:

            cxn_actual_spectrogram = cxn_reconstructed
            cxn_truth_spectrogram = truth_calculation
        else:
            cxn_actual_spectrogram = np.column_stack((cxn_actual_spectrogram, cxn_reconstructed))
            cxn_truth_spectrogram = np.column_stack(
                (cxn_truth_spectrogram, truth_calculation)
            )

    cxn_actual_axes.set_ylabel("Energy (GeV)")
    cxn_actual_axes.set_xlabel('Time (s)')
    cxn_actual_axes.set_title("CXN Reconstructed")
    __X_cxn, __Y_cxn = np.meshgrid(time_bins_x_axis/u.s, l_data[0]['Energy'])
    cxn_actual_axes.set_xlim(0.0001,20)
    cxn_actual_axes.set_ylim(0.01,0.1)

    cxn_actual_pc = cxn_actual_axes.pcolormesh(__X_cxn, __Y_cxn, cxn_actual_spectrogram)
    cxn_comp_fig.colorbar(cxn_actual_pc, ax=cxn_actual_axes, shrink=0.75, label=r'${cm}^{-2}$',
                            format='%.0e')
    cxn_actual_axes.set_xscale('log')
    # cxn_actual_axes.set_yscale('log')

    # cxn_truth_fig, cxn_truth_axes = plt.subplots(1,1,figsize=(8,8))
    cxn_truth_axes.set_ylabel('Energy (GeV)')
    cxn_truth_axes.set_xlabel('Time (s)')
    cxn_truth_axes.set_title('CXN Truth (SNOwGLoBES)')
    cxn_truth_axes.set_xlim(0.0001,20)
    cxn_truth_axes.set_ylim(0,0.1)
    __X_truth, __Y_truth = np.meshgrid(time_bins_x_axis/u.s, 10**np.array(ibd_cxn_actual['energy']))
    cxn_truth_pc = cxn_truth_axes.pcolormesh(__X_truth, __Y_truth, cxn_truth_spectrogram)
    cxn_comp_fig.colorbar(cxn_truth_pc, ax=cxn_truth_axes, shrink=0.75, label=r'${cm}^{-2}$',
                          format='%.0e')
    # cxn_truth_pc.set_clim(vmin=0, vmax=0.1)
    cxn_truth_axes.set_xscale('log')
    cxn_comp_fig.savefig('./time_cuts/cxn reconstruction comparison.png')
    cxn_comp_fig.show()

    one_slice_cxn_comp_fig, one_slice_cxn_truth_vs_recon = plt.subplots(1, 1, figsize=(8, 8))
    slice_index = 198

    one_slice_cxn_truth_vs_recon.scatter(
        l_data[0]['Energy'],
        np.transpose(cxn_actual_spectrogram)[slice_index],
        label="Reconstructed"
    )

    one_slice_cxn_truth_vs_recon.scatter(
        l_data[0]['Energy'],
        np.interp(
            l_data[0]['Energy'],
            10**np.array(ibd_cxn_actual['energy']),
            np.transpose(cxn_truth_spectrogram)[slice_index]),
        label="Truth"
    )

    one_slice_cxn_truth_vs_recon.set_xscale('log')
    one_slice_cxn_truth_vs_recon.set_yscale('log')

    # set titles
    one_slice_cxn_truth_vs_recon.set_title('Reconstructed vs Truth CXN in wc100kt30prct')

    # set axes titles
    one_slice_cxn_truth_vs_recon.set_xlabel('Energy (GeV)')

    one_slice_cxn_truth_vs_recon.set_ylabel(r'${cm}^{2}$')
    one_slice_cxn_comp_fig.legend()
    one_slice_cxn_comp_fig.savefig('./time_cuts/reconstruction over truth flux.png')
    one_slice_cxn_comp_fig.show()
