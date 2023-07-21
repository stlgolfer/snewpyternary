# file for striping unfolded files together into one
import click
import pandas as pd
import numpy as np
from meta_analysis import t_normalize
import snewpyternary as t

@click.command()
@click.option('-nux', required=True, help='Location of nux csv')
@click.option('-nue', required=True, help='Location of nue csv')
@click.option('-anue', required=True, help='Location of anue csv')
def bind(nux, nue, anue):
    nux_df = pd.read_csv(nux)
    nue_df = pd.read_csv(nue)
    anue_df = pd.read_csv(anue)

    # now that everything is loaded, need to put df into raw tuples
    raw_combined = list(zip(nux_df['unfolded'], anue_df['unfolded'], nue_df['unfolded']))
    ternary_points = t_normalize(raw_combined)
    print(ternary_points)
    t.create_default_flux_plot(ternary_points, "yuh", show=True)

if __name__ == '__main__':
    bind()
