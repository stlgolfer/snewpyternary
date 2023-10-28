# file for striping unfolded files together into one
import click
import pandas as pd
import numpy as np
import ternary

from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct
import snewpyternary as t
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict
import matplotlib.pyplot as plt
import math
import warnings

def project_ternary_points(points):
    scale = 100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    tax.set_title('Projection')
    # data is organized in top, right, left
    tax.bottom_axis_label(r'$\nu_x$')
    tax.right_axis_label(r'$\bar{\nu_e}$')
    tax.left_axis_label(r'$\nu_e$')
    tax.scatter(points)

    projected_points = tax.get_axes().collections[0].get_offsets().data
    # print(f'end: {projected_points[-1]}, intial: {projected_points[0]}')
    tax.show()
    # figure.close()
    return projected_points


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
    unfolded_pre_csum = list(
        zip(
            (nux_df['unfolded']/6),
            (anue_df['unfolded']),
            (nue_df['unfolded'])
        )
    )

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
            np.divide(nux_df['Ndet'], nux_df['dt']) / ebin,
            np.divide(anue_df['Ndet'], nux_df['dt']) / ebin,
            np.divide(nue_df['Ndet'], nux_df['dt']) / ebin
        )
    )

    ternary_points = t_normalize(unfolded_csum)
    unfolded_ternary_points_pre_csum = t_normalize(unfolded_pre_csum)

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

    #region time domain ternary error bars
    no_csum_frac_error_bars_td_nux = np.zeros(len(ternary_points))
    no_csum_frac_error_bars_td_anue = np.zeros(len(ternary_points))
    no_csum_frac_error_bars_td_nue = np.zeros(len(ternary_points))
    for i in range(len(no_csum_frac_error_bars_td_nux)):
        x = nux_df['unfolded'][i]/6
        y = (anue_df['unfolded'])[i]
        z = (nue_df['unfolded'])[i]
        A = list(zip(*ndet_raw_combined_per_time_pre_csum))[0][i]
        B = list(zip(*ndet_raw_combined_per_time_pre_csum))[1][i]
        C = list(zip(*ndet_raw_combined_per_time_pre_csum))[2][i]

        error_x = 3*math.sqrt((x**2*(B*C*(y + z)**2 + A*(C*y**2 + B*z**2)))/(A*B*C*(x + y + z)**4))
        error_y = 3*math.sqrt((y**2*(B*C*x**2 + A*(B*z**2 + C*(x + z)**2)))/(A*B*C*(x + y + z)**4))
        error_z = 3*math.sqrt(((B*C*x**2 + A*(C*y**2 + B*(x + y)**2))*z**2)/(A*B*C*(x + y + z)**4))

        no_csum_frac_error_bars_td_nux[i] = error_x*100
        no_csum_frac_error_bars_td_anue[i] = error_y*100
        no_csum_frac_error_bars_td_nue[i] = error_z*100
    #endregion

    #region time domain fractional error bars cumulative sum
    csum_frac_error_bars_td_nux = np.zeros(len(ternary_points))
    csum_frac_error_bars_td_anue = np.zeros(len(ternary_points))
    csum_frac_error_bars_td_nue = np.zeros(len(ternary_points))

    for i in range(len(csum_frac_error_bars_td_nux)):
        x = unfolded_transpose_csum[0][i]/6
        y = unfolded_transpose_csum[1][i]
        z = unfolded_transpose_csum[2][i]

        A = ndet_raw_transpose_csum[0][i]
        B = ndet_raw_transpose_csum[1][i]
        C = ndet_raw_transpose_csum[2][i]

        error_x = 3*math.sqrt((x**2*(B*C*(y + z)**2 + A*(C*y**2 + B*z**2)))/(A*B*C*(x + y + z)**4))
        error_y = 3*math.sqrt((y**2*(B*C*x**2 + A*(B*z**2 + C*(x + z)**2)))/(A*B*C*(x + y + z)**4))
        error_z = 3*math.sqrt(((B*C*x**2 + A*(C*y**2 + B*(x + y)**2))*z**2)/(A*B*C*(x + y + z)**4))

        csum_frac_error_bars_td_nux[i] = error_x*100
        csum_frac_error_bars_td_anue[i] = error_y*100
        csum_frac_error_bars_td_nue[i] = error_z*100
    #endregion

    # need to source error from original

    # print(ternary_points)
    # get the heatmap of it as well

    #region general plot settings
    ELINE_WIDTH = 1
    ECAP_SIZE = 4
    #endregion

    #region make some time domain plots as well
    time_domain_fig, (td_c_ax, td_c_frac_ax) = plt.subplots(1,2, figsize=(8,5))
    td_c_ax.errorbar(nux_df['time'], np.cumsum(nux_df['unfolded']/6), yerr=nux_csum_error_td, fmt='.', label=r'$\nu_x$ Unf.')
    td_c_ax.errorbar(nux_df['time'], np.cumsum(anue_df['unfolded']), yerr=anue_csum_error_td, fmt='.', label=r'$\bar{\nu_e}$ Unf.')
    td_c_ax.errorbar(nux_df['time'], np.cumsum(nue_df['unfolded']), yerr=nue_csum_error_td, fmt='.', label=r'$\nu_e$ Unf.')
    td_c_ax.set_xscale('log')
    td_c_ax.set_yscale('log')
    td_c_ax.set_xlabel('Mid-Point Time (s)')
    td_c_ax.set_ylabel(r'$\frac{neutrinos}{0.2*MeV*dt}$ Cumu.')
    td_c_ax.legend()

    td_c_ax_inset = td_c_ax.inset_axes([0.55,0.1, 0.4,0.5])
    td_c_ax_inset.errorbar(nux_df['time'], np.cumsum(nux_df['unfolded'] / 6), yerr=nux_csum_error_td, fmt='.',
                     label=r'$\nu_x$ Unf.')
    td_c_ax_inset.errorbar(nux_df['time'], np.cumsum(anue_df['unfolded']), yerr=anue_csum_error_td, fmt='.',
                     label=r'$\bar{\nu_e}$ Unf.')
    td_c_ax_inset.errorbar(nux_df['time'], np.cumsum(nue_df['unfolded']), yerr=nue_csum_error_td, fmt='.',
                     label=r'$\nu_e$ Unf.')
    td_c_ax_inset.set_xscale('log')

    td_c_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*ternary_points))[0],
        label=r'$\nu_x$ Unf.',
        yerr=csum_frac_error_bars_td_nux,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.'
    )
    td_c_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*ternary_points))[1],
        label=r'$\bar{\nu_e}$ Unf.',
        yerr=csum_frac_error_bars_td_anue,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.'
    )
    td_c_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*ternary_points))[2],
        label=r'$\nu_e$ Unf.',
        yerr=csum_frac_error_bars_td_nue,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.'
    )
    td_c_frac_ax.set_xscale('log')
    td_c_frac_ax.set_xlabel('Mid-Point Time (s)')
    td_c_frac_ax.set_ylabel('%')
    td_c_frac_ax.set_ylim(0, 100)
    td_c_frac_ax.legend()

    time_domain_fig.suptitle(title)
    time_domain_fig.savefig(f'./plots/{title} Cumulative Time Domain.png')
    time_domain_fig.show()

    time_domain_fig_no_csum, (td_ax, td_frac_ax) = plt.subplots(1,2,figsize=(8,5))

    print(no_csum_frac_error_bars_td_nux)
    td_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*unfolded_ternary_points_pre_csum))[0],
        yerr=no_csum_frac_error_bars_td_nux,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.',
        label=r'$\nu_x$ Unf.'
    )
    td_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*unfolded_ternary_points_pre_csum))[1],
        yerr=no_csum_frac_error_bars_td_anue,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.',
        label=r'$\bar{\nu_e}$ Unf.'
    )
    td_frac_ax.errorbar(
        nux_df['time'],
        list(zip(*unfolded_ternary_points_pre_csum))[2],
        yerr=no_csum_frac_error_bars_td_nue,
        elinewidth=ELINE_WIDTH,
        capsize=ECAP_SIZE,
        fmt='.',
        label=r'$\nu_e$ Unf.'
    )
    td_frac_ax.set_xscale('log')
    td_frac_ax.set_xlabel('Mid-Point Time (s)')
    td_frac_ax.set_ylim(0,100)
    td_frac_ax.set_ylabel('%')
    td_frac_ax.legend()

    td_ax.errorbar(nux_df['time'], nux_df['unfolded'] / 6, yerr=nux_error_td_pre_csum, fmt='.', label=r'$\nu_x$ Unf.')
    td_ax.errorbar(nux_df['time'], anue_df['unfolded'], yerr=anue_error_td_pre_csum, fmt='.', label=r'$\bar{\nu_e}$ Unf.')
    td_ax.errorbar(nux_df['time'], nue_df['unfolded'], yerr=nue_error_td_pre_csum, fmt='.', label=r'$\nu_e$ Unf.')
    td_ax.set_xscale('log')
    td_ax.set_yscale('log')
    td_ax.set_xlabel('Mid-Point Time (s)')
    td_ax.set_ylabel(r'$\frac{neutrinos}{0.2*MeV*dt}$')
    td_ax.legend()
    # let's also make the same plot but in regular y-scale
    td_ax_inset = td_ax.inset_axes([0.55,0.1, 0.4,0.5])
    td_ax_inset.errorbar(nux_df['time'], nux_df['unfolded'] / 6, yerr=nux_error_td_pre_csum, fmt='.', label=r'$\nu_x$ Unf.')
    td_ax_inset.errorbar(nux_df['time'], anue_df['unfolded'], yerr=anue_error_td_pre_csum, fmt='.',
                   label=r'$\bar{\nu_e}$ Unf.')
    td_ax_inset.errorbar(nux_df['time'], nue_df['unfolded'], yerr=nue_error_td_pre_csum, fmt='.', label=r'$\nu_e$ Unf.')
    td_ax_inset.set_xscale('log')

    time_domain_fig_no_csum.suptitle(title)
    time_domain_fig_no_csum.savefig(f'./plots/{title} Time Domain.png')
    time_domain_fig_no_csum.show()

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

    if heatmap:
        print("Generating heatmap (this might take a while)...")
        tax.heatmap(generate_heatmap_dict_phi_est(unfolded_csum, ternary_points, ndet_raw_combined_per_time, sigma_mult=3),
                    cmap=plt.get_cmap('PiYG'))
        print("Done")
    tax.plot_colored_trajectory(t_normalize(unfolded_csum), cmap=plt.get_cmap('binary'))

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

    #region calculate curliness of ternary cumul plot
    # ternary points is this
    # max_point = None
    # max_distance = 0
    # main_line = ternary_distance(ternary_points[0], ternary_points[-1])
    # for point in ternary_points[1:-1]:
    #     a = ternary_distance(point, ternary_points[-1])
    #     b = ternary_distance(point, ternary_points[0])
    #     s = (a + main_line + b) / 2  # Heron's semi-perimeter
    #     h = math.sqrt(4 * s * (s - a) * (s - b) * (s - main_line) / main_line ** 2)
    #     if h > max_distance:
    #         max_distance = h
    #         max_point = point
    # main_line_slope = (ternary_points[0][1] - ternary_points[-1][1]) / (ternary_points[0][0] - ternary_points[-1][0])
    #
    # main_line_eqn_output = main_line_slope * (max_point[0] - ternary_points[0][0]) - ternary_points[0][1]  # point slope
    # if main_line_slope < 0 and max_point[1] < main_line_eqn_output:
    #     max_distance = max_distance * -1
    # elif main_line_slope > 0 and max_point[1] > main_line_eqn_output:
    #     max_distance = max_distance * -1
    # # now max distance is the "curl" with sign correction
    # print(f'Curliness is: {max_distance}')

    # let's try getting the projected points directly from matplotlib

    # now we also want to write this out to a file. we will just keep appending to the same file since binder is
    # typically used repeatedly
    # calculate slope using projected data
    proj_points = project_ternary_points(ternary_points)
    proj_f = proj_points[-1]
    proj_i = proj_points[0]
    slope = (proj_f[1]-proj_i[1])/(proj_f[0]-proj_i[0])
    print(f'slope is {slope}')
    curliness_file = open("curliness_from_binder.csv", "a")
    curliness_file.write(f'{title},{slope}\n')
    curliness_file.close()
    #endregion

if __name__ == '__main__':
    bind()
