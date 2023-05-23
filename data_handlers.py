#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 23:16:27 2022
This contains data handlers for handling individual model data
@author: phyics
"""

import numpy as np
from abc import ABC, abstractmethod

supported_detectors = [
    'scint20kt',
    'ar40kt',
    'wc100kt30prct'
    ]

__NA__: float = 6.0221E23

class __DetectorProxyConfiguration__(ABC):
    '''
    This is an abstract class that can be used to describe which channels should be used as proxies for neutrino
    flavors.
    '''


    @abstractmethod
    def h_scint20kt(self) -> ([str], [str], [str]):
        '''
        This is the abstract handler for the scint20kt
        Parameters
        ----------
        data: a dictionary with channels as keys

        Returns
        -------
        The summed event rates for the time bin
        '''
        pass

    @abstractmethod
    def h_ar40kt(self) -> ([str], [str], [str]):
        '''
                This is the abstract handler for the scint20kt
                Parameters
                ----------
                data: a dictionary with channels as keys

                Returns
                -------
                The summed event rates for the time bin
                '''
        pass

    @abstractmethod
    def h_wc100kt30prct(self) -> ([str], [str], [str]):
        '''

        Returns
        -------
        For each proxy, there is a list of channels that are to be summed for calculating the proxy
        '''
        pass

    def Nt_scint20kt(self) -> [float]:
        '''
        Nt (number of targets calculation used in a simple unfolding)

        Returns
        -------
        arr: array of Nt values in order of nux, nue, a_nue
        '''
        return [1.0, 1.0, 1.0]

    def Nt_ar40kt(self) -> [float]:
        '''
        Nt (number of targets calculation used in a simple unfolding)

        Returns
        -------
        arr: array of Nt values in order of nux, nue, a_nue
        '''
        return [1.0, 1.0, 1.0]

    def Nt_wc100kt30prct(self):
        '''
        Nt (number of targets calculation used in a simple unfolding)

        Returns
        -------
        arr: array of Nt values in order of nux, nue, a_nue
        '''

        return [1.0, 1.0, 1.0]

    @abstractmethod
    def __str__(self):
        '''

        Returns
        -------
        the name of the proxy configuration
        '''
        pass

    '''
    Data handler axes labels definitions
    '''

    def same_axes(self):
        return [r'$\nu_x$ Proxy', r'$\nu_e$ Proxy', r'$\bar{\nu_e}$ Proxy']

    def axes_scint20kt(self):
        return self.same_axes()

    def axes_ar40kt(self):
        return self.same_axes()

    def axes_wc100kt30prct(self):
        return self.same_axes()

    def build_detector_profiles(self):
        '''
        Builds a profile for a detector that includes a dictionary of data handlers
        and the corresponding data labels (used for labeling axes)

        Returns
        -------
        detector_defs : dictionary
            A dictionary with the supported_detectors as the keys and the respective values are the handler and axes
            defs
            TODO: might not have to run this at runtime: might be able to just construct it
        '''

        detector_defs = {}
        module_scope = self # __import__(__name__)
        # unfortunately, now the handler code just got more complicated, because now we have to generically sum
        # channels based on the names
        for detector in supported_detectors:
            # print(f'h_{detector}')
            h_channels = getattr(module_scope, f'h_{detector}')()
            # print(h_channels)

            detector_defs[detector] = {
                'chans_to_add': h_channels,
                'axes': getattr(module_scope, f'axes_{detector}'),
                'N_t': getattr(module_scope, f'Nt_{detector}')
            }
        return detector_defs

class ConfigAggregateDetectors(__DetectorProxyConfiguration__):
    '''
        What follows are the calculations that are made for each time bin when collated
        data comes through. Each time bin has event rates per energy bin. Order should
        be nc, nu_e+es, then inverse beta decay events.

        Treat all the NC as nux,  all the nue-nucleus and ES as nue, and IBD+nuebar-nucleus as nuebar.

        Each function here should return nux, nue, and anue, respectively if using the same_axes label
        '''

    def h_scint20kt(self):
        # must return a list of a, b, c
        # nue_plus_es=np.sum(data['nue_C12'])+np.sum(data['nue_C13']+data['e']) # nue
        # ibd = np.sum(data['ibd'])
        # nc = np.sum(data['nc'])
        return (['nc'], ['nue_C12', 'nue_C13', 'e'], ['ibd'])

    def h_ar40kt(self):
        # nue_plus_es = np.sum(data['nue_Ar40'])+np.sum(data['e'])
        # nc = np.sum(data['nc'])
        # nue_bar = np.sum(data['nuebar_Ar40'])
        return (['nc'], ['nue_Ar40', 'e'], ['nuebar_Ar40'])

    def h_wc100kt30prct(self):
        # ibd = np.sum(data['ibd'])
        # nue_plus_es=np.sum(data['nue_O16'])+np.sum(data['e'])
        # nc = np.sum(data['nc'])
        return (['nc'], ['nue_O16', 'e'], ['ibd'])

    def __str__(self):
        return "AgDet"

class ConfigBestChannel(__DetectorProxyConfiguration__):

    def h_scint20kt(self):
        # ibd = np.sum(data['ibd'])
        # nc = np.sum(data['nc'])
        # return [nc, 0, ibd]
        return (['nc'], [], ['ibd'])

    def h_ar40kt(self):
        return ([], ['nue_Ar40', 'e'], [])

    def h_wc100kt30prct(self):
        return ([], [], ['ibd'])

    def __str__(self):
        return "BstChnl"

    def Nt_scint20kt(self) -> [float]:
        calc = 50*(100E9)*__NA__/12
        return [calc, 1, 1]

    def Nt_argon40kt(self) -> [float]:
        return [1, 40*(100E9)*__NA__/39.9, 1]

    def Nt_wc100kt30prct(self):
        return [1, 1, (100)*(100E9)*2*__NA__/9]
