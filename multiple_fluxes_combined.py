import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

nakazato_vanilla = pd.read_csv('./flux_saves/Nakazato_2013_s0_NoTransformation_BstChnl_flux_save_log.csv')
nakazato_imo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_flux_save.csv')
nakazato_nmo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_NMO_BstChnl_flux_save.csv')

unfolded_nakazato_vanilla = pd.read_csv('./binders/Nakazato_2013_s0_NoTransformation_nocsum_unfolded.csv')
unfolded_nakazato_imo = pd.read_csv('./binders/Nakazato_2013_s0_AdiabaticMSW_IMO_nocsum_unfolded.csv')
unfolded_nakazato_nmo = pd.read_csv('./binders/Nakazato_2013_s0_AdiabaticMSW_NMO_nocsum_unfolded.csv')
fig, [(ax1, ax2, ax3), (ax4, ax5, ax6)] = plt.subplots(2, 3, figsize=(11,9))

times = np.array(nakazato_vanilla['time'])
times = times + np.abs(times[0]) + 0.001
# print(times)

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

ax1.plot(times, nakazato_vanilla['raw_data_nux']/4, label=r'$\nu_x$')
ax1.plot(times, nakazato_vanilla['raw_data_nue'], label=r'$\nu_e$')
ax1.plot(times, nakazato_vanilla['raw_data_anue'], label=r'$\bar{\nu}_e$')
ax1.set_xscale('log')
# ax1.set_xlim(0.001,1)
ax1.set_ylabel(r'Neutrinos/($cm^2$*Time Bin)', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax1.set_title('NoTrans. Flux')
ax1.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax1.legend(fontsize=14, loc='upper left')
ax1_bins = slice_bins(ax1, times)
# test cumsum
ax1_csum = ax1.twinx()
ax1_csum.plot(
    times,
    np.cumsum(nakazato_vanilla['raw_data_nux']/4),
    linestyle='dashed',
    alpha=0.5
)
ax1_csum.plot(
    times,
    np.cumsum(nakazato_vanilla['raw_data_nue']),
    linestyle='dashed',
    alpha=0.5
)
ax1_csum.plot(
    times,
    np.cumsum(nakazato_vanilla['raw_data_anue']),
    linestyle='dashed',
    alpha=0.5
)

log_times = np.array(nakazato_imo['time'])
log_times = log_times + np.abs(log_times[0]) + 0.001
# print(log_times)

ax2.plot(log_times, nakazato_imo['raw_data_nux']/4, label=r'$\nu_x$')
ax2.plot(log_times, nakazato_imo['raw_data_nue'], label=r'$\nu_e$')
ax2.plot(log_times, nakazato_imo['raw_data_anue'], label=r'$\bar{\nu}_e$')
ax2.set_xscale('log')
# ax2.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax2.set_title('Adia. IMO Flux')
ax2.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax2.legend(fontsize=14)
slice_bins(ax2, log_times)

ax2_csum = ax2.twinx()
ax2_csum.plot(
    times,
    np.cumsum(nakazato_imo['raw_data_nux']/4),
    linestyle='dashed',
    alpha=0.5
)
ax2_csum.plot(
    times,
    np.cumsum(nakazato_imo['raw_data_nue']),
    linestyle='dashed',
    alpha=0.5
)
ax2_csum.plot(
    times,
    np.cumsum(nakazato_imo['raw_data_anue']),
    linestyle='dashed',
    alpha=0.5
)


ax3.plot(log_times, nakazato_nmo['raw_data_nux']/4, label=r'$\nu_x$')
ax3.plot(log_times, nakazato_nmo['raw_data_nue'], label=r'$\nu_e$')
ax3.plot(log_times, nakazato_nmo['raw_data_anue'], label=r'$\bar{\nu}_e$')
# td_ax_inset = td_ax.inset_axes([0.55,0.1, 0.4,0.5])
ax3.set_title('Adia. NMO Flux')
ax3.set_xscale('log')
# ax3.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax3.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax3.legend(fontsize=14)
slice_bins(ax3, log_times)

ax3_csum = ax3.twinx()
ax3_csum.plot(
    times,
    np.cumsum(nakazato_nmo['raw_data_nux']/4),
    linestyle='dashed',
    alpha=0.5
)
ax3_csum.plot(
    times,
    np.cumsum(nakazato_nmo['raw_data_nue']),
    linestyle='dashed',
    alpha=0.5
)
ax3_csum.plot(
    times,
    np.cumsum(nakazato_nmo['raw_data_anue']),
    linestyle='dashed',
    alpha=0.5
)
ax3_csum.set_ylabel('Cumulative Sum', fontsize=12)
print(f'Showing every {len(times)/ax1_bins} bins, each of which are {(times[-1]-times[0])/len(times)} s')


# now create the unfolded versions
ax4.plot(times, unfolded_nakazato_vanilla['nux_df']*(3/2), label=r'$\nu_x$')
ax4.plot(times, unfolded_nakazato_vanilla['nue_df'], label=r'$\nu_e$')
ax4.plot(times, unfolded_nakazato_vanilla['anue_df'], label=r'$\bar{\nu_e}$')
ax4.set_xscale('log')
# ax4.set_xlim(0.001,1)
ax4.set_ylabel(r'Neutrinos/($cm^2$*Time Bin)', fontsize=16)
# ax4.set_title('Neutronization Burst', fontsize=14)
ax4.set_title('NoTran. Unfolded')
ax4.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax4.legend(fontsize=14)
ax4_bins = slice_bins(ax4, times)

ax4_csum = ax4.twinx()
ax4_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_vanilla['nux_df']*(3/2)),
    linestyle='dashed',
    alpha=0.5
)
ax4_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_vanilla['nue_df']),
    linestyle='dashed',
    alpha=0.5
)
ax4_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_vanilla['anue_df']),
    linestyle='dashed',
    alpha=0.5
)

log_times = np.array(unfolded_nakazato_imo['time'])
log_times = log_times + np.abs(log_times[0]) + 0.001
# print(log_times)

ax5.plot(log_times, unfolded_nakazato_imo['nux_df']*(3/2), label=r'$\nu_x$')
ax5.plot(log_times, unfolded_nakazato_imo['nue_df'], label=r'$\nu_e$')
ax5.plot(log_times, unfolded_nakazato_imo['anue_df'], label=r'$\bar{\nu_e}$')
ax5.set_xscale('log')
# ax2.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax5.set_title('Adia. IMO Unfolded')
ax5.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax5.legend(fontsize=14)
slice_bins(ax5, log_times)

ax5_csum = ax5.twinx()
ax5_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_imo['nux_df']*(3/2)),
    linestyle='dashed',
    alpha=0.5
)
ax5_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_imo['nue_df']),
    linestyle='dashed',
    alpha=0.5
)
ax5_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_imo['anue_df']),
    linestyle='dashed',
    alpha=0.5
)



ax6.plot(log_times, unfolded_nakazato_nmo['nux_df']*(3/2), label=r'$\nu_x$')
ax6.plot(log_times, unfolded_nakazato_nmo['nue_df'], label=r'$\nu_e$')
ax6.plot(log_times, unfolded_nakazato_nmo['anue_df'], label=r'$\bar{\nu_e}$')
# td_ax_inset = td_ax.inset_axes([0.55,0.1, 0.4,0.5])
ax6.set_title('Adia. NMO Unfolded')
ax6.set_xscale('log')
# ax3.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax6.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax6.legend(fontsize=14)
slice_bins(ax6, log_times)

ax6_csum = ax6.twinx()
ax6_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_nmo['nux_df']*(3/2)),
    linestyle='dashed',
    alpha=0.5
)
ax6_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_nmo['nue_df']),
    linestyle='dashed',
    alpha=0.5
)
ax6_csum.plot(
    times,
    np.cumsum(unfolded_nakazato_nmo['anue_df']),
    linestyle='dashed',
    alpha=0.5
)
ax6_csum.set_ylabel('Cumulative Sum', fontsize=12)

fig.tight_layout()
fig.savefig('./flux_comparison_all.png')
