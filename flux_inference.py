import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
# from math import exp
from scipy.special import gamma

if __name__ == '__main__':
    def garching(x, n0, alpha, ev):
        apo = alpha + 1 #"alpha plus one"
        # print(type(n0))
        alphas = np.power(apo,apo)
        denom = ev * gamma(apo)
        frac = x/ev
        return n0*(alphas/denom)*(np.power(frac, alpha))*np.exp(-apo*frac)

    fig, ax = plt.subplots(1, 1)
    df = pd.read_csv('./flux_spectra/Nakazato_2013_s0_AdiabaticMSW_IMO_BstChnl_nu_e_bar_flux_spectra.csv')
    energy_bins = np.array(df['flux_energy_bins']*1000)
    flux_vs_energy = np.array(df['flux_vs_energy'])

    # now attempt to fit using the garching parameters
    (N0, Alpha0, Ev), pcov = curve_fit(
        garching,
        energy_bins,
        flux_vs_energy,
        p0=(np.max(flux_vs_energy), 1e-3, 50),
        bounds=(
            (0,0, 0),#lower
            (np.Inf, 10, 1000)#upper?
        )
        # full_output=True
        # lowerBounds = (0, -np.Inf, 0)
    )
    print(Alpha0)
    
    ax.plot(energy_bins, garching(energy_bins, N0, Alpha0, Ev))
    ax.scatter(energy_bins, flux_vs_energy)
    ax.set_ylabel(r'$\frac{dn}{dE_\nu}$')
    ax.set_xlabel(r'$E_\nu$' + ' (MeV)')
    ax.set_yscale('log')
    fig.savefig('./flux_spectra.png')