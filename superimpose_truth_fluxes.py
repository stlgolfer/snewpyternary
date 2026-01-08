import click
import pandas as pd
import numpy as np
import ternary

from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct, process_flux
import sys
sys.path.insert(0,'./SURF2020fork')
from SURF2020fork.ternary_helpers import generate_heatmap_dict_phi_est, generate_heatmap_dict
import matplotlib.pyplot as plt
import math
import warnings

import data_handlers
import snewpyternary as t
from model_wrappers import snewpy_models, sn_model_default_time_step

def plot_flux(tax, ordering, model_number=0, cmap='binary'):
    config = t.MetaAnalysisConfig(
        snewpy_models['Nakazato_2013'],
        [model_number],
        f'AdiabaticMSW_{ordering}',
        proxy_config=data_handlers.ConfigBestChannel()
    )

    _, raw_data, _ = process_flux(config, model_number)
    raw_data_processed = [(x[0] / 6, x[1], x[2]) for x in raw_data]
    # flux_td_fig, flux_td_tax = t.create_default_flux_plot(t_normalize(np.cumsum(np.array(raw_data_processed),0)))
    plotting_data = t_normalize(np.cumsum(np.array(raw_data_processed), 0))

    tax.plot_colored_trajectory(plotting_data, cmap=plt.get_cmap(cmap))
    tax.scatter([tuple(plotting_data[-1])], marker='s', color='yellow', linewidth=10)
    tax.scatter([tuple(plotting_data[0])], marker='^', linewidth=10, color='cyan')
    # tax.savefig('./test flux.png')
    # flux_td_tax.show()

if __name__ == '__main__':
    TERNARY_AXES_LABEL_FONT_SIZE = 20
    scale = 100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    # tax.set_title(rf'Ndet')
    # data is organized in top, right, left
    tax.bottom_axis_label(r'$\nu_x$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.right_axis_label(r'$\bar{\nu_e}$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.left_axis_label(r'$\nu_e$', fontsize=TERNARY_AXES_LABEL_FONT_SIZE)
    tax.ticks(axis='lbr', linewidth=1, multiple=scale / 10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    # paint_track(tax, 'NMO', 'copper')
    # paint_track(tax, 'IMO', 'binary')
    plot_flux(tax, 'IMO', model_number=0, cmap='binary')
    plot_flux(tax, 'NMO', model_number=0, cmap='copper')
    plot_flux(tax, 'IMO', model_number=2, cmap='binary')
    plot_flux(tax, 'NMO', model_number=2, cmap='copper')

    tax.savefig('Nakazato_2013 truth fluxes superimposed.png')