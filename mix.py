from tkinter import filedialog as fd
import tkinter as tk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from meta_analysis import t_normalize
from snewpyternary import clean_newline
import ternary
import click
import glob
from tqdm import tqdm

def process_mixture(cxn_filename, zeta_filename, show=True):
    # now that we have filenames, we have to load the data
    cxn = pd.read_csv(cxn_filename)
    zeta = pd.read_csv(zeta_filename)

    # plot the cxn in linear time
    cxn_plot, (cxn_nux_ax, cxn_anue_ax, cxn_nue_ax) = plt.subplots(1,3)
    cxn_nux_ax.plot(cxn['time_bins'], cxn[r'$\nu_x$ Proxy'])
    cxn_nux_ax.set_title('nux Proxy')

    cxn_anue_ax.plot(cxn['time_bins'], cxn[r'$\bar{\nu_e}$ Proxy'])
    cxn_anue_ax.set_title('anue Proxy')

    cxn_nue_ax.plot(cxn['time_bins'], cxn[r'$\nu_e$ Proxy'])
    cxn_nue_ax.set_title('nue Proxy')
    cxn_plot.show()

    # want to sample from the cxn as much as possible using the zeta bin number coordinates
    cxn_bin_numbers = [x for x in range(len(cxn['time_bins']))]
    zeta_bin_numbers = [x for x in range(len(zeta['time_bins']))]

    cxn_nux_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\nu_x$ Proxy'])
    cxn_anue_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\bar{\nu_e}$ Proxy'])
    cxn_nue_resampled = np.interp(zeta_bin_numbers, cxn_bin_numbers, cxn[r'$\nu_e$ Proxy'])

    nux_phi_est = np.divide(zeta[r'$\nu_x$ Proxy'], cxn_nux_resampled)
    anue_phi_est = np.divide(zeta[r'$\bar{\nu_e}$ Proxy'], cxn_anue_resampled)
    nue_phi_est = np.divide(zeta[r'$\nu_e$ Proxy'], cxn_nue_resampled)

    # now make the plot
    phi_est_fig, (phi_nux_axes, phi_anue_axes, phi_nue_axes) = plt.subplots(1, 3)
    phi_nux_axes.plot(zeta_bin_numbers, nux_phi_est, label='nux est.')
    phi_anue_axes.plot(zeta_bin_numbers, anue_phi_est, label='anue est.')
    phi_nue_axes.plot(zeta_bin_numbers, nue_phi_est, label='nue est.')
    phi_est_fig.legend()
    if show:
        phi_est_fig.show()

    # and plot in ternary space?
    phi_zipped = list(zip(nux_phi_est, anue_phi_est, nue_phi_est))
    phi_t_normalized = t_normalize(phi_zipped)

    scale = 100
    figure, tax = ternary.figure(scale=scale)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=scale / 10)
    # tax.set_title(title)
    # data is organized in top, right, left
    # apparently this is in flux formatting
    tax.bottom_axis_label('nux')
    tax.right_axis_label('nuebar')
    tax.left_axis_label('nue')

    widths = np.linspace(0.01, 1, num=len(phi_t_normalized))
    colorid = 0
    for p in range(len(phi_t_normalized) - 1):
        if (p + 1 >= len(phi_t_normalized)):
            break
        tax.line(phi_t_normalized[p], phi_t_normalized[p + 1], color=(
            widths[p] if colorid == 0 else 0, widths[p] if colorid == 1 else 0, widths[p] if colorid == 2 else 0, 1),
                 linestyle=':', linewidth=3)
    # tax.scatter(phi_t_normalized)

    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')  # disables regular matlab plot axes
    # print(phi_t_normalized)
    tax_title = f'CXN: {cxn_filename.split("/")[-1]}\n Zeta: {zeta_filename.split("/")[-1]}'
    tax.set_title(tax_title)
    tax.savefig(f'./mixtures/{clean_newline(tax_title)}.png')
    if show:
        tax.show()

@click.command('Mix')
@click.argument('mode', required=True, type=str, nargs=1)
def start(mode):
    if mode=='permute':
        print("Permute mode selected")
        # first load and list all cxns
        cxn_files = glob.glob('./cxns/*.csv')
        zeta_files = glob.glob('./zetas/*.csv')
        with tqdm(total=len(cxn_files)*len(zeta_files)) as pbar:
            for cxn in cxn_files:
                for zeta in zeta_files:
                    process_mixture(cxn, zeta, show=False)
                    pbar.update(1)


    elif mode=='select':
        print("Select mode")
        root = tk.Tk()
        root.withdraw()
        cxn_filename = fd.askopenfilename(title="Open Cross-Section File", initialdir='./cxns')
        print(f'Selected CXN: {cxn_filename}')
        zeta_filename = fd.askopenfilename(title="Open Zeta File", initialdir='./zetas')
        print(f'Selected Zeta: {zeta_filename}')
        process_mixture(cxn_filename, zeta_filename)
    else:
        print("Hard code mode")
        cxn_filename = "/home/acm123/snewpy-snowglobes-ternary/cxns/CXNs for Nakazato_2013 AdiabaticMSW_IMO BstChnlRED nakazato-shen-z0.004-t_rev100ms-s20.0.fits Logged Bins .csv"
        zeta_filename = "/home/acm123/snewpy-snowglobes-ternary/zetas/zetas for Nakazato_2013 AdiabaticMSW_IMO BstChnlRED nakazato-shen-z0.004-t_rev100ms-s20.0.fits Logged Bins .csv"
        process_mixture(cxn_filename, zeta_filename)



if __name__ == '__main__':
    print(f"Modes are:\n=============================\n 1) select: pick the files \n 2) permute: permute ./cxns and ./zetas")
    start()
