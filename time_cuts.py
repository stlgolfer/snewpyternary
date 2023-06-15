import snewpyternary as t
import snowglobes_wrapper
from model_wrappers import snewpy_models, sn_model_default_time_step
import data_handlers
from meta_analysis import process_flux, process_detector
import matplotlib.pyplot as plt
import numpy as np

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
    plot_data, raw_data, l_data, zeta, sigma_average = process_detector(config,0,'wc100kt30prct')

    # want 0.05s, 0.07s, 0.1s, 0.3s, 0.5s, 1, 3, 5, 7, 10
    time_cut_indexes = [114,117,121,135,142,153,170,178,184,190]
    for time in time_cut_indexes:
        fig, (flux_axes, ndet_axes) = plt.subplots(1,2,figsize=(16,8))
        # want to superimpose flux and Ndet, so make sure GeV
        flux_energy_spectra = np.linspace(0, 100, 501)  # * MeV  # 1MeV
        flux_axes.scatter(flux_energy_spectra,labeled[time][4],label='Flux')
        flux_axes.set_title('Flux')
        flux_axes.set_xlabel('Energy (MeV)')
        flux_axes.set_ylabel('neutrinos/(cm^2 * MeV)')

        ndet_axes.scatter(l_data[time]['Energy']*1000,l_data[time]['ibd'], label='Ndet')
        ndet_axes.set_title('Ndet')
        ndet_axes.set_xlabel('Energy (MeV)')
        ndet_axes.set_ylabel('Event count / MeV')

        fig.suptitle(rf'Flux vs Ndet Spect. $t\approx{time_bins_x_axis[time]}$')

        fig.savefig(f'./time_cuts/{time_bins_x_axis[time]}.png')
