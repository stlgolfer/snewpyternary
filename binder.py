# file for striping unfolded files together into one
import click
import pandas as pd
import numpy as np
import ternary

from meta_analysis import t_normalize
import snewpyternary as t
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict
import matplotlib.pyplot as plt
import warnings

@click.command()
@click.option('-nux', required=True, help='Location of nux csv')
@click.option('-nue', required=True, help='Location of nue csv')
@click.option('-anue', required=True, help='Location of anue csv')
@click.option('--title', required=False, help='Title of the ternary plot', default='')
@click.option('--heatmap', required=False, help='Create heatmap?', default=True)
def bind(nux, nue, anue, title, heatmap):
    nux_df = pd.read_csv(nux)
    nue_df = pd.read_csv(nue)
    anue_df = pd.read_csv(anue)
    # time = nux_df['time']

    # now that everything is loaded, need to put df into raw tuples
    unfolded_csum = list(
        zip(
            np.cumsum(nux_df['unfolded']/6),
            np.cumsum(anue_df['unfolded']),
            np.cumsum(nue_df['unfolded'])
        )
    )
    unfolded_pre_csum = t_normalize(list(
        zip(
            (nux_df['unfolded']/6),
            (anue_df['unfolded']),
            (nue_df['unfolded'])
        )
    ))

    ebin = 0.02e-3
    ndet_raw_combined_per_time = list(
        zip(
            np.cumsum(np.divide(nux_df['Ndet'],nux_df['dt'])/ebin),
            np.cumsum(np.divide(anue_df['Ndet'],nux_df['dt'])/ebin),
            np.cumsum(np.divide(nue_df['Ndet'],nux_df['dt'])/ebin)
        )
    )

    ndet_raw_combined_per_time_pre_csum = list(
        zip(
            (np.divide(nux_df['Ndet'], nux_df['dt']) / ebin),
            (np.divide(anue_df['Ndet'], nux_df['dt']) / ebin),
            (np.divide(nue_df['Ndet'], nux_df['dt']) / ebin)
        )
    )

    #region time domain error bars
    unfolded_transpose_csum = list(zip(*unfolded_csum))
    ndet_raw_transpose_csum = list(zip(*ndet_raw_combined_per_time))

    nux_csum_error_td = 3*np.divide(unfolded_transpose_csum[0], np.sqrt(ndet_raw_transpose_csum[0]))
    anue_csum_error_td = 3*np.divide(unfolded_transpose_csum[1], np.sqrt(ndet_raw_transpose_csum[1]))
    nue_csum_error_td = 3*np.divide(unfolded_transpose_csum[2], np.sqrt(ndet_raw_transpose_csum[2]))

    unfolded_transpose_pre_csum = list(zip(*unfolded_pre_csum))
    ndet_raw_transpose_pre_csum = list(zip(*ndet_raw_combined_per_time_pre_csum))

    nux_error_td_pre_csum = 3*np.divide(unfolded_transpose_pre_csum[0], np.sqrt(ndet_raw_transpose_pre_csum[0]))
    anue_error_td_pre_csum = 3*np.divide(unfolded_transpose_pre_csum[1], np.sqrt(ndet_raw_transpose_pre_csum[1]))
    nue_error_td_pre_csum = 3*np.divide(unfolded_transpose_pre_csum[2], np.sqrt(ndet_raw_transpose_pre_csum[2]))
    #endregion

    # need to source error from original
    ternary_points = t_normalize(unfolded_csum)
    # print(ternary_points)
    # get the heatmap of it as well

    #region make some time domain plots as well
    time_domain_fig, ((td_c_ax, td_c_frac_ax), (td_ax, td_frac_ax)) = plt.subplots(2,2, figsize=(8,8))
    td_c_ax.errorbar(nux_df['time'], np.cumsum(nux_df['unfolded']/6), yerr=nux_csum_error_td, fmt='o', label=r'$\nu_x$ Unf.')
    td_c_ax.errorbar(nux_df['time'], np.cumsum(anue_df['unfolded']), yerr=anue_csum_error_td, fmt='o', label=r'$\bar{\nu_e}$ Unf.')
    td_c_ax.errorbar(nux_df['time'], np.cumsum(nue_df['unfolded']), yerr=nue_csum_error_td, fmt='o', label=r'$\nu_e$ Unf.')
    td_c_ax.set_xscale('log')
    td_c_ax.set_yscale('log')
    td_c_ax.set_xlabel('Mid-Point Time (s)')
    td_c_ax.set_ylabel(r'$\frac{neutrinos}{0.2*MeV*dt}$ Cumu.')
    td_c_ax.legend()

    td_c_frac_ax.scatter(nux_df['time'], list(zip(*ternary_points))[0], label=r'$\nu_x Unf.')
    td_c_frac_ax.scatter(nux_df['time'], list(zip(*ternary_points))[1], label=r'$\bar{\nu_e} Unf.')
    td_c_frac_ax.scatter(nux_df['time'], list(zip(*ternary_points))[2], label=r'$\nu_e Unf.')
    td_c_frac_ax.set_xscale('log')
    td_c_frac_ax.set_xlabel('Mid-Point Time (s)')
    td_c_frac_ax.set_ylabel('%')
    td_c_frac_ax.legend()

    td_frac_ax.scatter(nux_df['time'], list(zip(*unfolded_pre_csum))[0], label=r'$\nu_x Unf.')
    td_frac_ax.scatter(nux_df['time'], list(zip(*unfolded_pre_csum))[1], label=r'$\bar{\nu_e} Unf.')
    td_frac_ax.scatter(nux_df['time'], list(zip(*unfolded_pre_csum))[2], label=r'$\nu_e Unf.')
    td_frac_ax.set_xscale('log')
    td_frac_ax.set_xlabel('Mid-Point Time (s)')
    td_frac_ax.set_ylabel('%')
    td_frac_ax.legend()

    td_ax.errorbar(nux_df['time'], nux_df['unfolded'] / 6, yerr=nux_error_td_pre_csum, fmt='o', label=r'$\nu_x$ Unf.')
    td_ax.errorbar(nux_df['time'], anue_df['unfolded'], yerr=anue_error_td_pre_csum, fmt='o', label=r'$\bar{\nu_e}$ Unf.')
    td_ax.errorbar(nux_df['time'], nue_df['unfolded'], yerr=nue_error_td_pre_csum, fmt='o', label=r'$\nu_e$ Unf.')
    td_ax.set_xscale('log')
    td_ax.set_yscale('log')
    td_ax.set_xlabel('Mid-Point Time (s)')
    td_ax.set_ylabel(r'$\frac{neutrinos}{0.2*MeV*dt}$')
    td_ax.legend()
    time_domain_fig.suptitle(title)
    time_domain_fig.savefig(f'./plots/{title} Time Domain.png')
    time_domain_fig.show()

    #endregion

    #region ternary diagram for phi_est_flux
    scale = 100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    tax.set_title(rf'{title} $\phi_e$')
    # data is organized in top, right, left
    tax.bottom_axis_label(r'$\nu_x$')
    tax.right_axis_label(r'$\bar{\nu_e}$')
    tax.left_axis_label(r'$\nu_e$')

    # tax.scatter(points=plotting_data, color="red")
    widths = np.linspace(0.01, 1, num=len(unfolded_csum))
    # for p in range(len(plotting_data) - 1):
    #     if (p + 1 >= len(plotting_data)):
    #         break
        # tax.line(plotting_data[p], plotting_data[p + 1], color=(widths[p], 0, 0, 1), linestyle=':', linewidth=3)
    if heatmap:
        print("Generating heatmap (this might take a while)...")
        tax.heatmap(generate_heatmap_dict_phi_est(unfolded_csum, ternary_points, ndet_raw_combined_per_time, sigma_mult=3),
                    cmap=plt.get_cmap('PiYG'))
        print("Done")
    colormap = {}
    for i,p in enumerate(t_normalize(unfolded_csum)):
        colormap[tuple(p)] = widths[i]
    tax.plot_colored_trajectory(t_normalize(unfolded_csum), cmap=plt.get_cmap('binary'))
    # tax.heatmap(colormap, cmap=plt.get_cmap('Greys'))
    # hm, what if we tried stacking heatmaps?


    tax.ticks(axis='lbr', linewidth=1, multiple=scale / 10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')

    # fig, tax = t.create_default_flux_plot(t_normalize(unfolded_csum), rf'{title} $\phi_e$', save=False, show=False)
    tax.show()
    tax.savefig(f'./fluxes/{title} Unfolded.png')
    #endregion

    Ndet_fig, Ndet_tax = t.create_default_flux_plot(t_normalize(ndet_raw_combined_per_time), f'{title} Ndet',save=False,show=False)
    if heatmap:
        Ndet_tax.heatmap(generate_heatmap_dict(ndet_raw_combined_per_time,t_normalize(ndet_raw_combined_per_time)))
    Ndet_tax.show()
    Ndet_tax.savefig(f'./plots/{title} Ndet.png')

if __name__ == '__main__':
    bind()
