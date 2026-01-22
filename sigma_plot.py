import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

if __name__ == '__main__':
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 5))

    def paint_channel(ax, fname, chan_name):
        # read in an example. there will be a cross section for each channel
        df = pd.read_csv(fname)
        sigma_avg = df['sigma average']
        time = df['time']
        dt = df['dt']
        ax.plot(time, sigma_avg, label='<cxn>. (previous)')
        ax.set_xscale('log')
        # now also plot the avg version of this
        avg = (np.array(sigma_avg) @ np.array(dt)) / np.sum(dt)
        print(avg)
        ax.plot(time, avg*np.ones_like(time), label='Avg.')
        ax.set_title(chan_name)
        ax.legend()
        ax.set_xlabel(r'$ t+t_0 $')
    
    ax1.set_ylabel(r'$ \sim\sigma~{cm}^2 $')# a bit unsure of the units
    fig.suptitle('Nakazato s0 <cxn> Comparison')
    paint_channel(ax1, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_e_bar_sigma_average.csv', r'$\bar{\nu}_e$' + ' channel')
    paint_channel(ax2, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_e_sigma_average.csv', r'$\nu_e$' + ' channel')
    paint_channel(ax3, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_mu_sigma_average.csv', r'$\nu_\mu$' + ' channel')
    fig.savefig('cxn.png')