import matplotlib.pyplot as plt
import numpy as np

from meta_analysis import process_detector
# use_log=False

import data_handlers
import snewpyternary as t
from model_wrappers import snewpy_models, sn_model_default_time_step

if __name__ == '__main__':
    
    config = t.MetaAnalysisConfig(
        snewpy_models['Nakazato_2013'],
        [0],
        'NoTransformation',
        proxy_config=data_handlers.ConfigBestChannel()
    )
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4.8))
    argon_pc = process_detector(config, 0, 'ar40kt', ax=ax1)
    water_pc = process_detector(config, 0, 'wc100kt30prct', ax=ax2)
    scint_pc = process_detector(config, 0, 'scint20kt', ax=ax3)

    fig.colorbar(argon_pc, ax=ax1, shrink=0.75, location='right', format='%.0e')
    fig.colorbar(water_pc, ax=ax2, shrink=0.75, location='right', format='%.0e')
    fig.colorbar(scint_pc, ax=ax3, shrink=0.75, location='right', format='%.0e').set_label(fontsize=14,label='Event Count / (0.2 MeV * Time Bin)')

    ax1.set_xlabel(r'Time + $t_0$ (s)', fontsize=13)
    ax2.set_xlabel(r'Time + $t_0$ (s)', fontsize=13)
    ax3.set_xlabel(r'Time + $t_0$ (s)', fontsize=13)

    ax1.set_ylabel('Energy (MeV)', fontsize=14)

    ax1.set_title('Argon 40kt', fontsize=13)
    ax2.set_title('Water 100kt', fontsize=13)
    ax3.set_title('Scintillator 20kt', fontsize=13)

    ax1.tick_params(axis='x', labelsize=11)
    ax2.tick_params(axis='x', labelsize=11)
    ax3.tick_params(axis='x', labelsize=11)

    ax1.tick_params(axis='y', labelsize=11)
    ax2.tick_params(axis='y', labelsize=11)
    ax3.tick_params(axis='y', labelsize=11)

    fig.tight_layout()
    fig.savefig(f'./spectra_signal_figure.png')