import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def slice_bins(axes, xvals, slices=20, color='black', alpha=0.1):
    # want to plot vertical lines where bins are
    axes_y_lims = axes.get_ylim()
    # print(len(xvals))
    bar_slices = np.floor(np.linspace(0, len(xvals) - 1, slices)).astype(int)
    # print(bar_slices)
    axes.vlines(
        xvals[bar_slices],
        axes_y_lims[0], axes_y_lims[1],
        linestyles='dashed',
        color=color,
        alpha=alpha
        )
    print(f'Showing every {len(xvals)/len(bar_slices)} bins, each of which are {(xvals[-1]-xvals[0])/len(xvals)} s')
    return len(bar_slices)

if __name__ == '__main__':
    nakazato_vanilla = pd.read_csv('./flux_saves/Nakazato_2013_s0_NoTransformation_BstChnl_flux_save_log.csv')
    nakazato_imo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_flux_save.csv')
    nakazato_nmo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_NMO_BstChnl_flux_save.csv')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,6))

    times = np.array(nakazato_vanilla['time'])
    times = times + np.abs(times[0]) + 0.001
    # print(times)

    ax1.plot(times, nakazato_vanilla['raw_data_nux']/4, label=r'$\nu_x$')
    ax1.plot(times, nakazato_vanilla['raw_data_nue'], label=r'$\nu_e$')
    ax1.plot(times, nakazato_vanilla['raw_data_anue'], label=r'$\bar{\nu}_e$')
    ax1.set_xscale('log')
    # ax1.set_xlim(0.001,1)
    ax1.set_ylabel(r'Neutrinos/($cm^2$*Time Bin)', fontsize=16)
    # ax1.set_title('Neutronization Burst', fontsize=14)
    ax1.set_title('NoTran. Log Time Bins')
    ax1.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
    ax1.legend(fontsize=14)
    ax1_bins = slice_bins(ax1, times)

    log_times = np.array(nakazato_imo['time'])
    log_times = log_times + np.abs(log_times[0]) + 0.001
    # print(log_times)

    ax2.plot(log_times, nakazato_imo['raw_data_nux']/4, label=r'$\nu_x$')
    ax2.plot(log_times, nakazato_imo['raw_data_nue'], label=r'$\nu_e$')
    ax2.plot(log_times, nakazato_imo['raw_data_anue'], label=r'$\bar{\nu}_e$')
    ax2.set_xscale('log')
    # ax2.set_ylabel(r'neutrinos/cm^2', fontsize=16)
    # ax1.set_title('Neutronization Burst', fontsize=14)
    ax2.set_title('Adia. IMO Log Time Bins')
    ax2.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
    ax2.legend(fontsize=14)
    slice_bins(ax2, log_times)


    ax3.plot(log_times, nakazato_nmo['raw_data_nux']/4, label=r'$\nu_x$')
    ax3.plot(log_times, nakazato_nmo['raw_data_nue'], label=r'$\nu_e$')
    ax3.plot(log_times, nakazato_nmo['raw_data_anue'], label=r'$\bar{\nu}_e$')
    # td_ax_inset = td_ax.inset_axes([0.55,0.1, 0.4,0.5])
    ax3.set_title('Adia. NMO Log Time Bins')
    ax3.set_xscale('log')
    # ax3.set_ylabel(r'neutrinos/cm^2', fontsize=16)
    # ax1.set_title('Neutronization Burst', fontsize=14)
    ax3.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
    ax3.legend(fontsize=14)
    slice_bins(ax3, log_times)
    print(f'Showing every {len(times)/ax1_bins} bins, each of which are {(times[-1]-times[0])/len(times)} s')
    fig.tight_layout()
    fig.savefig('./flux_comparison.png')
