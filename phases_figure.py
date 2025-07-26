import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
ax1.set_ylabel(r'neutrinos/cm^2', fontsize=16)
ax1.set_title('Neutronization Burst', fontsize=14)
ax1.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax1.legend(fontsize=14)

ax2.plot(times, df['raw_data_nux']/6, label=r'$\nu_x$')
ax2.plot(times, df['raw_data_nue'], label=r'$\nu_e$')
ax2.plot(times, df['raw_data_anue'], label=r'$\bar{\nu_e}$')
# ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.8, 20)
ax2.set_ylim(1e7,1e9)
ax2.set_title('Cooling', fontsize=14)
ax2.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax2.legend(fontsize=14)
fig.savefig('./phases.png')