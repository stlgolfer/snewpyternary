from tkinter import filedialog as fd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from meta_analysis import t_normalize
import ternary

ALLOW_GUI = True

cxn_filename = "/home/acm123/snewpy-snowglobes-ternary/cxns/CXNs for Nakazato_2013 AdiabaticMSW_IMO BstChnlRED nakazato-shen-z0.004-t_rev100ms-s20.0.fits Logged Bins .csv"
zeta_filename = "/home/acm123/snewpy-snowglobes-ternary/zetas/zetas for Nakazato_2013 AdiabaticMSW_IMO BstChnlRED nakazato-shen-z0.004-t_rev100ms-s20.0.fits Logged Bins .csv"

if ALLOW_GUI:
    cxn_filename = fd.askopenfilename(title="Open Cross-Section File", initialdir='./cxns')
    print(f'Selected CXN: {cxn_filename}')
    zeta_filename = fd.askopenfilename(title="Open Zeta File", initialdir='./zetas')
    print(f'Selected Zeta: {zeta_filename}')

# now that we have filenames, we have to load the data
cxn = pd.read_csv(cxn_filename)
zeta = pd.read_csv(zeta_filename)

# want to sample from the cxn as much as possible using the zeta bin number coordinates
cxn_bin_numbers = [x for x in range(len(cxn['time_bins']))]
zeta_bin_numbers = [x for x in range(len(zeta['time_bins']))]

cxn_nux_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\nu_x$ Proxy'])
cxn_anue_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\bar{\nu_e}$ Proxy'])
cxn_nue_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\nu_e$ Proxy'])

nux_phi_est = np.divide(zeta[r'$\nu_x$ Proxy'],cxn_nux_resampled)
anue_phi_est = np.divide(zeta[r'$\bar{\nu_e}$ Proxy'],cxn_anue_resampled)
nue_phi_est = np.divide(zeta[r'$\nu_e$ Proxy'],cxn_nue_resampled)

# now make the plot
phi_est_fig, (phi_nux_axes, phi_anue_axes, phi_nue_axes) = plt.subplots(1,3)
phi_nux_axes.plot(zeta_bin_numbers, nux_phi_est, label='nux est.')
phi_anue_axes.plot(zeta_bin_numbers, anue_phi_est, label='anue est.')
phi_nue_axes.plot(zeta_bin_numbers, nue_phi_est, label='nue est.')
phi_est_fig.legend()
phi_est_fig.show()

# and plot in ternary space?
phi_zipped = list(zip(nux_phi_est, anue_phi_est, nue_phi_est))
phi_t_normalized = t_normalize(phi_zipped)

scale=100
figure, tax = ternary.figure(scale=scale)
tax.boundary(linewidth=2.0)
tax.gridlines(color="blue", multiple=scale/10)
# tax.set_title(title)
# data is organized in top, right, left
# apparently this is in flux formatting
tax.bottom_axis_label('nux')
tax.right_axis_label('nuebar')
tax.left_axis_label('nue')
tax.scatter(phi_t_normalized)

tax.clear_matplotlib_ticks()
tax.get_axes().axis('off') # disables regular matlab plot axes
print(phi_t_normalized)
tax.set_title(f'CXN: {cxn_filename.split("/")[-1]}\nZeta: {zeta_filename.split("/")[-1]}')
tax.show()
