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
def bind(nux, nue, anue, title):
    nux_df = pd.read_csv(nux)
    nue_df = pd.read_csv(nue)
    anue_df = pd.read_csv(anue)
    # time = nux_df['time']

    # now that everything is loaded, need to put df into raw tuples
    raw_combined = list(
        zip(
            np.cumsum(nux_df['unfolded']/6),
            np.cumsum(anue_df['unfolded']),
            np.cumsum(nue_df['unfolded'])
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

    # need to source error from original
    ternary_points = t_normalize(raw_combined)
    # print(ternary_points)
    # get the heatmap of it as well

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
    widths = np.linspace(0.01, 1, num=len(raw_combined))
    # for p in range(len(plotting_data) - 1):
    #     if (p + 1 >= len(plotting_data)):
    #         break
        # tax.line(plotting_data[p], plotting_data[p + 1], color=(widths[p], 0, 0, 1), linestyle=':', linewidth=3)
    print("Generating heatmap (this might take a while)...")
    tax.heatmap(generate_heatmap_dict_phi_est(raw_combined, ternary_points, ndet_raw_combined_per_time, sigma_mult=3),
                cmap=plt.get_cmap('PiYG'))
    print("Done")
    colormap = {}
    for i,p in enumerate(t_normalize(raw_combined)):
        colormap[tuple(p)] = widths[i]
    tax.plot_colored_trajectory(t_normalize(raw_combined), cmap=plt.get_cmap('binary'))
    # tax.heatmap(colormap, cmap=plt.get_cmap('Greys'))
    # hm, what if we tried stacking heatmaps?


    tax.ticks(axis='lbr', linewidth=1, multiple=scale / 10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')

    # fig, tax = t.create_default_flux_plot(t_normalize(raw_combined), rf'{title} $\phi_e$', save=False, show=False)
    tax.show()
    tax.savefig(f'./fluxes/{title} Unfolded.png')

    Ndet_fig, Ndet_tax = t.create_default_flux_plot(t_normalize(ndet_raw_combined_per_time), f'{title} Ndet',save=False,show=False)
    Ndet_tax.heatmap(generate_heatmap_dict(ndet_raw_combined_per_time,t_normalize(ndet_raw_combined_per_time)))
    Ndet_tax.show()
    Ndet_tax.savefig(f'./plots/{title} Ndet.png')

if __name__ == '__main__':
    bind()
