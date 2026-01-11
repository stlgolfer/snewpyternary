import click
import pandas as pd
import numpy as np
import ternary

from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict
import matplotlib.pyplot as plt
import math
import warnings

import data_handlers
import snewpyternary as t
from model_wrappers import snewpy_models, sn_model_default_time_step

# @click.command()
# @click.option('-nux', required=True, help='Location of nux csv')
# @click.option('-nue', required=True, help='Location of nue csv')
# @click.option('-anue', required=True, help='Location of anue csv')
# @click.option('--title', required=False, help='Title of the ternary plot', default='')
# @click.option('--heatmap', required=False, help='Create heatmap?', default=True)
TERNARY_AXES_LABEL_FONT_SIZE = 20
def bind(tax, nux, nue, anue, cmap=None, heat=True):
    nux_df = pd.read_csv(nux)
    nue_df = pd.read_csv(nue)
    anue_df = pd.read_csv(anue)
    ndet_raw_combined_per_time = list(
        zip(
            np.cumsum(nux_df['Ndet']),
            np.cumsum(anue_df['Ndet']),
            np.cumsum(nue_df['Ndet'])
        )
    )
    ternary_points = t_normalize(ndet_raw_combined_per_time)
    if heat:
        tax.heatmap(
            generate_heatmap_dict(
                ndet_raw_combined_per_time,
                ternary_points
            ),
            cmap=plt.get_cmap('PuOr'),
            colorbar=False
        )
    tax.plot_colored_trajectory(ternary_points, cmap=plt.get_cmap(cmap))
    tax.scatter([tuple(ternary_points[-1])], marker='s', color='yellow', linewidth=10)
    tax.scatter([tuple(ternary_points[0])], marker='^', linewidth=10, color='cyan')
    
    # Ndet_fig, Ndet_tax = t.create_default_flux_plot(t_normalize(ndet_raw_combined_per_time_pre_csum), f'Ndet',save=False,show=False)

    # Ndet_tax.show()
    # Ndet_tax.savefig(f'./test Ndet.png')

def paint_track(tax, i, ordering, cmap=None, heat=False):
    model = 'Nakazato_2013'
    bind(tax,
        f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_mu_sigma_average.csv',
        f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_e_sigma_average.csv',
        f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_e_bar_sigma_average.csv',
        cmap=cmap,
        heat=heat
        )
    tax.set_title(ordering, fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    
def set_diagram(tax):
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    # tax.set_title(rf'Ndet')
    # data is organized in top, right, left
    tax.bottom_axis_label('Neutral Current', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.right_axis_label('IBD', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.left_axis_label('Charged Current', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.ticks(axis='lbr', linewidth=1, multiple=scale / 10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    tax._redraw_labels()

if __name__ == '__main__':
    # python binder.py -nux ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_AdiabaticMSW_NMO --heatmap=${heat}
    # we want to make a 1x2 grid that has an NMO and IMO to show that the error regions really are quite bad
    # we can make another one in the appendix?
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    with_heatmap = True

    scale = 100
    fig1, tax1 = ternary.figure(ax=ax1, scale=scale)
    set_diagram(tax1)
    paint_track(tax1, 0, 'NMO', 'binary', heat=with_heatmap)

    fig2, tax2 = ternary.figure(ax=ax2, scale=scale)
    set_diagram(tax2)
    paint_track(tax2, 0, 'IMO', 'binary', heat=with_heatmap)
    # fig1.savefig('test ndet.png')
    # fig1.savefig('test ndet.png') # is apparently needed to render the figures
    # fig2.savefig('test ndet.png')
    fig.savefig('Nakazato_2013 ndet superimposed.png')
