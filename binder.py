# file for striping unfolded files together into one
import click
import pandas as pd
import numpy as np
from meta_analysis import t_normalize
import snewpyternary as t
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict
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
            np.cumsum(nux_df['unfolded']),
            np.cumsum(anue_df['unfolded']),
            np.cumsum(nue_df['unfolded'])
        )
    )
    ndet_raw_combined = list(zip(nux_df['Ndet'], anue_df['Ndet'], nue_df['Ndet']))

    # need to source error from original
    ternary_points = t_normalize(raw_combined)
    # print(ternary_points)
    # get the heatmap of it as well

    fig, tax = t.create_default_flux_plot(ternary_points, title, save=False, show=False)
    print("Generating heatmap (this might take a while)...")
    tax.heatmap(generate_heatmap_dict(ndet_raw_combined, ternary_points), cmap=plt.get_cmap('PiYG'))
    print("Done")
    tax.show()

if __name__ == '__main__':
    bind()
