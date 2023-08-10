# file for striping unfolded files together into one
import click
import pandas as pd
import numpy as np
from meta_analysis import t_normalize
import snewpyternary as t
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict, generate_heatmap_dict_v2
import matplotlib.pyplot as plt

@click.command()
@click.option('-nux', required=True, help='Location of nux csv')
@click.option('-nue', required=True, help='Location of nue csv')
@click.option('-anue', required=True, help='Location of anue csv')
@click.option('--title', required=False, help='Title of the ternary plot', default='')
def bind(nux, nue, anue, title):
    nux_df = pd.read_csv(nux)
    nue_df = pd.read_csv(nue)
    anue_df = pd.read_csv(anue)

    # now that everything is loaded, need to put df into raw tuples
    raw_combined = list(
        zip(
            nux_df['unfolded'],
            anue_df['unfolded'],
            nue_df['unfolded']
        )
    )

    # delta_time = np.subtract(nux_df['time'][1:], nux_df['time'][0:-1])
    ndet_raw_combined = list(
        zip(
            nux_df['Ndet'],
            anue_df['Ndet'],
            nue_df['Ndet']
        )
    )

    ebin = 0.02e-3
    ndet_raw_combined_per_time = list(
        zip(
            np.divide(nux_df['Ndet'],nux_df['dt'])/ebin,
            np.divide(anue_df['Ndet'],nux_df['dt'])/ebin,
            np.divide(nue_df['Ndet'],nux_df['dt'])/ebin
        )
    )

    # need to source error from original
    ternary_points = t_normalize(raw_combined)
    # print(ternary_points)
    # get the heatmap of it as well

    fig, tax = t.create_default_flux_plot(ternary_points, title, save=False, show=False)
    print("Generating heatmap (this might take a while)...")
    tax.heatmap(generate_heatmap_dict_phi_est(raw_combined, ternary_points, ndet_raw_combined_per_time), cmap=plt.get_cmap('PiYG'))
    print("Done")
    tax.show()

    Ndet_fig, Ndet_tax = t.create_default_flux_plot(t_normalize(ndet_raw_combined_per_time), "",save=False,show=False)
    Ndet_tax.heatmap(generate_heatmap_dict(ndet_raw_combined_per_time,t_normalize(ndet_raw_combined_per_time)))
    # ndet_heatmap = generate_heatmap_dict_v2(ndet_raw_combined,t_normalize(ndet_raw_combined))
    # Ndet_tax.heatmap(ndet_heatmap, cmap=plt.get_cmap('PiYG'))
    Ndet_tax.show()

if __name__ == '__main__':
    bind()
