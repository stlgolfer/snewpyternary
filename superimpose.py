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

# change matplotlib font sizes
TERNARY_AXES_LABEL_FONT_SIZE = 20

@click.command()
@click.option('--config', required=True, help='File path to superimpose configuration')
@click.option('--title', required=False, help='Title of the ternary plot', default='')
def superimpose(config, title):
    # plan here will be to start a diagram
    # then load the list of bound files and tax on the figure

    #region ternary diagram for phi_est_flux
    scale = 100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    tax.set_title(rf'{title} $\phi_e$')
    # data is organized in top, right, left
    tax.bottom_axis_label(r'$\nu_x$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.right_axis_label(r'$\bar{\nu_e}$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.left_axis_label(r'$\nu_e$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)

    with open(config, 'r') as file:
        files = file.readlines()
        for f in files:
            df = pd.read_csv(f)
            unfolded_csum = list(
                zip(
                    df['nux_df'],
                    df['anue_df'],
                    df['nue_df']
                )
            )
            ternary_points = t_normalize(unfolded_csum)
            tax.plot_colored_trajectory(t_normalize(unfolded_csum), cmap=plt.get_cmap('binary'))


            # we also want to paint the starting point and end point differently
            tax.scatter([tuple(ternary_points[-1])], marker='s', color='yellow', linewidth=10)
            tax.scatter([tuple(ternary_points[0])], marker='^', linewidth=10, color='cyan')

    tax.ticks(axis='lbr', linewidth=1, multiple=scale / 10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')

    # fig, tax = t.create_default_flux_plot(t_normalize(unfolded_csum), rf'{title} $\phi_e$', save=False, show=False)
    tax.show()
    tax.savefig(f'./fluxes/test.png')
    print(f'Ternary diagram painted {len(ternary_points)} points')
    #endregion

if __name__ == '__main__':
    superimpose()
