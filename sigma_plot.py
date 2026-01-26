import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import data_handlers

if __name__ == '__main__':
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 5))

    def paint_channel(ax, fname, chan_name, nt):
        # read in an example. there will be a cross section for each channel
        df = pd.read_csv(fname)
        sigma_avg = df['sigma average']
        time = df['time']
        dt = df['dt']
        ax.plot(time, sigma_avg, label='<cxn>. (previous)')
        ax.set_xscale('log')
        unfolded = df['unfolded']
        ndet = df['Ndet']

        # now also plot the avg version of this
        avg = (np.array(sigma_avg) @ np.array(dt)) / np.sum(dt)
        print(avg)
        ax.plot(time, avg*np.ones_like(time), label='Avg.')
        ax.plot(time, np.divide(ndet, unfolded)/nt, linestyle='dashed')
        ax.set_title(chan_name)
        ax.legend()
        ax.set_xlabel(r'$ t+t_0 $')
    
    ax1.set_ylabel(r'$ \sim\sigma~{cm}^2 $')# a bit unsure of the units
    fig.suptitle('Nakazato s0 <cxn> Comparison')
    nts = data_handlers.ConfigBestChannel()
    paint_channel(ax1, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_e_bar_sigma_average.csv', r'$\bar{\nu}_e$' + ' channel', nts.Nt_wc100kt30prct()[2])
    paint_channel(ax2, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_e_sigma_average.csv', r'$\nu_e$' + ' channel', nts.Nt_ar40kt()[1])
    paint_channel(ax3, './sigmas/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_mu_sigma_average.csv', r'$\nu_\mu$' + ' channel', nts.Nt_scint20kt()[0])
    fig.savefig('cxn.png')