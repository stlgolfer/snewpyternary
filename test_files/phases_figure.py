import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

if __name__ == '__main__':
    # fig, (burst_ax, accretion_ax, cooling_ax) = plt.subplots(1,3,figsize=(10,8))
    fig, burst_ax = plt.subplots(1,1)

    df = pd.read_csv('../flux_saves/Nakazato_2013_s0_NoTransformation_BstChnl_flux_save.csv')
    # need to correct time axis
    time_raw = np.array(df['time'])
    offset = 0.00001
    time = np.abs(time_raw[0]) + time_raw+offset
    print(f'Time offset is {time_raw[0]}')

    #TODO: needed to divide by 6 for nux
    burst_ax.scatter(time, np.array(df['raw_data_nux'])/4,label=r'$\nu_x$') # needs to be divide by four here since we're using truth flux
    burst_ax.scatter(time, df['raw_data_anue'], label=r'$\bar{\nu}_e$')
    burst_ax.scatter(time, df['raw_data_nue'], label=r'$\nu_e$')
    burst_ax.set_ylabel(r'$\frac{neutrinos}{cm^2}$')
    burst_ax.legend()
    burst_ax.set_xlabel(r'Time $+ t_0$')

    # for neutronization burst, log bins
    # burst_ax.set_title('Neutronization Burst')
    # burst_ax.set_xscale('log')
    # burst_ax.set_xlim(0.02, 0.35)

    # for accretion, linear bins
    # burst_ax.set_title('Accretion')
    # burst_ax.set_xlim(0.07, 2)

    # for cooling, linear bins
    burst_ax.set_xlim(14,20)
    burst_ax.set_ylim(0,0.5e9)
    burst_ax.set_title('Cooling')

    # accretion_ax.scatter(time, np.array(df['raw_data_nux'])/6, label=r'$\nu_x$')
    # accretion_ax.scatter(time, df['raw_data_anue'], label=r'$\bar{\nu}_e$')
    # accretion_ax.scatter(time, df['raw_data_nue'], label=r'$\nu_e$')
    # accretion_ax.set_xscale('log')
    # accretion_ax.set_ylabel(r'$\frac{neutrinos}{cm^2}$')
    # accretion_ax.set_title('Accretion')
    # accretion_ax.legend()
    # accretion_ax.set_xlabel(r'Time $+ t_0$')
    # accretion_ax.set_xlim(time_raw[0],1.2)
    #
    # cooling_ax.scatter(time, df['raw_data_nux'], label=r'$\nu_x$')
    # cooling_ax.scatter(time, df['raw_data_anue'], label=r'$\bar{\nu}_e$')
    # cooling_ax.scatter(time, df['raw_data_nue'], label=r'$\nu_e$')
    # cooling_ax.set_xscale('log')
    # cooling_ax.set_ylabel(r'$\frac{neutrinos}{cm^2}$')
    # cooling_ax.set_title('Cooling')
    # cooling_ax.legend()
    # cooling_ax.set_xlabel(r'Time $+ t_0$')
    # cooling_ax.set_xlim(10,20)
    # cooling_ax.set_ylim(0,1e10)
    fig.show()