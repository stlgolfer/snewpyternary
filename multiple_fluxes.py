import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

nakazato_vanilla = pd.read_csv('./flux_saves/Nakazato_2013_s0_NoTransformation_BstChnl_flux_save.csv')
nakazato_imo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_flux_save.csv')
nakazato_nmo = pd.read_csv('./flux_saves/Nakazato_2013_s0_AdiabaticMSW_NMO_BstChnl_flux_save.csv')
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,6))

times = np.array(nakazato_vanilla['time'])
times = times + np.abs(times[0]) + 0.001
print(times)

ax1.plot(times, nakazato_vanilla['raw_data_nux']/4, label=r'$\nu_x$')
ax1.plot(times, nakazato_vanilla['raw_data_nue'], label=r'$\nu_e$')
ax1.plot(times, nakazato_vanilla['raw_data_anue'], label=r'$\bar{\nu_e}$')
ax1.set_xscale('log')
# ax1.set_xlim(0.001,1)
ax1.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax1.set_title('NoTran. Uniform Time Bins')
ax1.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax1.legend(fontsize=14)

log_times = np.array(nakazato_imo['time'])
log_times = log_times + np.abs(log_times[0]) + 0.001
print(log_times)

ax2.plot(log_times, nakazato_imo['raw_data_nux']/4, label=r'$\nu_x$')
ax2.plot(log_times, nakazato_imo['raw_data_nue'], label=r'$\nu_e$')
ax2.plot(log_times, nakazato_imo['raw_data_anue'], label=r'$\bar{\nu_e}$')
ax2.set_xscale('log')
# ax2.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax2.set_title('Adia. IMO Log Time Bins')
ax2.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax2.legend(fontsize=14)

ax3.plot(log_times, nakazato_nmo['raw_data_nux']/4, label=r'$\nu_x$')
ax3.plot(log_times, nakazato_nmo['raw_data_nue'], label=r'$\nu_e$')
ax3.plot(log_times, nakazato_nmo['raw_data_anue'], label=r'$\bar{\nu_e}$')
# td_ax_inset = td_ax.inset_axes([0.55,0.1, 0.4,0.5])
ax3.set_title('Adia. NMO Log Time Bins')
ax3.set_xscale('log')
# ax3.set_ylabel(r'neutrinos/cm^2', fontsize=16)
# ax1.set_title('Neutronization Burst', fontsize=14)
ax3.set_xlabel(r'Time $t+t_0$ (s)', fontsize=12)
ax3.legend(fontsize=14)
fig.tight_layout()
fig.savefig('./flux_comparison.png')
