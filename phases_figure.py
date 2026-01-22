import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from multiple_fluxes import slice_bins

df = pd.read_csv('./flux_saves/Nakazato_2013_s0_NoTransformation_BstChnl_flux_save.csv')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4.8))
times = np.array(df['time'])

times = times + np.abs(times[0]) + 0.001
print(times)
ax1.plot(times, df['raw_data_nux']/4, label=r'$\nu_x$')
ax1.plot(times, df['raw_data_nue'], label=r'$\nu_e$')
ax1.plot(times, df['raw_data_anue'], label=r'$\bar{\nu_e}$')
ax1.set_xscale('log')
ax1.set_xlim(4e-2,1)
ax1.set_ylabel(r'Neutrinos/($cm^2$*Time Bin)', fontsize=16)
ax1.set_title('Neutronization Burst', fontsize=14)
ax1.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax1.legend(fontsize=14)
ax1_bins = slice_bins(ax1, times, slices=1000)
print(f'LEFT: showing every {len(times)/ax1_bins} bins, each of which are {(times[-1]-times[0])/len(times)} s')
print(f'So the spacing is {(times[-1]-times[0])/ax1_bins}')

ax2.plot(times, df['raw_data_nux']/4, label=r'$\nu_x$')
ax2.plot(times, df['raw_data_nue'], label=r'$\nu_e$')
ax2.plot(times, df['raw_data_anue'], label=r'$\bar{\nu_e}$')
# ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.8, 20)
ax2.set_ylim(1e7,1e9)
ax2.set_title('Cooling', fontsize=14)
ax2.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax2.legend(fontsize=14)
ax2_bins = slice_bins(ax2, times)
print(f'RIGHT: showing every {len(times)/ax2_bins} bins, each of which are {(times[-1]-times[0])/len(times)} s')
print(f'So the spacing is {(times[-1]-times[0])/ax2_bins}')
fig.savefig('./phases.png')