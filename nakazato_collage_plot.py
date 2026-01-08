import matplotlib.pyplot as plt
import numpy as np
from binder import TERNARY_AXES_LABEL_FONT_SIZE, WATER_NDET_MULTIPLIER
import ternary
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict
import pandas as pd
from meta_analysis import t_normalize

class BinderTriple:
    def __init__(self, model, i, ordering):
        nux_df = pd.read_csv(
            f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_mu_sigma_average.csv'
        )
        nue_df = pd.read_csv(
            f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_e_sigma_average.csv'
        )
        anue_df = pd.read_csv(
            f'./sigmas/{model}_s{i}_AdiabaticMSW_{ordering}_BstChnl_nu_e_bar_sigma_average.csv'
        )

        # now calculate as binder does
        self.unfolded_csum = list(
            zip(
                np.cumsum(nux_df['unfolded']/6),
                np.cumsum(anue_df['unfolded']),
                np.cumsum(nue_df['unfolded'])
            )
        )
        self.unfolded_pre_csum = list(
            zip(
                (nux_df['unfolded']/6),
                (anue_df['unfolded']),
                (nue_df['unfolded'])
            )
        )
        self.ndet_raw_combined_per_time = list(
        zip(
            np.cumsum(nux_df['Ndet']),
            np.cumsum(anue_df['Ndet']*WATER_NDET_MULTIPLIER),
            np.cumsum(nue_df['Ndet'])
        )
        )

        self.ndet_raw_combined_per_time_pre_csum = list(
            zip(
                nux_df['Ndet'],
                anue_df['Ndet']*WATER_NDET_MULTIPLIER,
                nue_df['Ndet']
            )
        )
        

def make_ternary_plot(
    ax,
    unfolded_csum,
    ndet_raw_combined_per_time,
    heatmap=True,
    title=None
    # unfolded_pre_csum,
    # unfolded_ternary_points_pre_csum,
    # ndet_raw_combined_per_time_pre_csum,
    ):

    scale = 100
    scale = 100
    figure, tax = ternary.figure(ax=ax, scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10, linewidth=1)
    # tax.set_title(rf'{title} $\phi_e$')
    if title is not None:
        tax.set_title(title)
    # data is organized in top, right, left
    tax.bottom_axis_label(r'$\nu_x$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.right_axis_label(r'$\bar{\nu_e}$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.left_axis_label(r'$\nu_e$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)

    
    print("Generating heatmap (this might take a while)...")
    ternary_points = t_normalize(unfolded_csum)
    if heatmap:
        tax.heatmap(generate_heatmap_dict_phi_est(unfolded_csum, ternary_points, ndet_raw_combined_per_time),
                cmap=plt.get_cmap('PuOr'), colorbar=False)
    print("Done")
    tax.plot_colored_trajectory(ternary_points, cmap=plt.get_cmap('binary'))
    # we also want to paint the starting point and end point differently
    tax.scatter([tuple(ternary_points[-1])], marker='s', color='yellow', linewidth=10)
    tax.scatter([tuple(ternary_points[0])], marker='^', linewidth=10, color='cyan')

    tax.ticks(axis='lbr', linewidth=1, multiple=100)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    # fig.savefig('test.png')
    # ax = tax

def make_four_plot(start, end):
    frac = 1

    fig = plt.figure(figsize=(frac*10,frac*11), constrained_layout=True) # 8, 15 seems to work
    ax = fig.subplot_mosaic(
        [
            [f's{start}IMO', f's{start}NMO'],
            [f's{end}IMO', f's{end}NMO'],
            # ['s3IMO', 's3NMO'],
            # ['s4IMO', 's4NMO']
        ],
        gridspec_kw={
            'width_ratios':[1,1],
            'height_ratios': [1,1,]
        },
        subplot_kw={
            'box_aspect':1
        }
    )
    blank = False
    for s in [start, end]:
        for o in ['IMO', 'NMO']:
            # test a triple
            if not blank:
                binder = BinderTriple(
                    'Nakazato_2013',
                    s,
                    o
                )
                make_ternary_plot(
                    ax[f's{s}{o}'],
                    binder.unfolded_csum,
                    binder.ndet_raw_combined_per_time,
                    heatmap=True,
                    # title=f's{s}'
                )
            else:
                continue
    # fig.constrained()
    fig.savefig(f'./mosaic{start}_{end}.png')

if __name__ == '__main__':
    make_four_plot(0,2)
    make_four_plot(3,4)