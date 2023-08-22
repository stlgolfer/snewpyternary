#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 23:43:12 2022

@author: physics
"""

from snewpy.models import Bollig_2016, Fornax_2019, Fornax_2021, Kuroda_2020, Nakazato_2013, Sukhbold_2015, Tamborra_2014, Walk_2018, Walk_2019, Warren_2020, Zha_2021
from astropy import units as u

class SNEWPYModel:
    """
    SNEWPY models have several available configs, so it's convenient to just
    bundle this information togeter
    """

    def __init__(self,file_paths: list, model_type: str, model: any):
        self.file_paths: list = file_paths
        self.model_type: str = model_type
        self.model: any = model

snewpy_models = {}

__sn_model_def_time_step_lib = {}
#TODO: might make a function that pulls the time and serves a default time if a key isn't found
__sn_model_def_time_step_lib['Bollig_2016'] = 0.1*u.s #0.014
__sn_model_def_time_step_lib['Fornax_2019'] = 0.1*u.s
__sn_model_def_time_step_lib['Fornax_2021'] = 0.025*u.s # don't use for now as of Issue #20
__sn_model_def_time_step_lib['Kuroda_2020'] = 0.005*u.s
__sn_model_def_time_step_lib['Nakazato_2013'] = 0.1*u.s # 0.1*u.s
__sn_model_def_time_step_lib['Sukhbold_2015'] = 0.1*u.s
__sn_model_def_time_step_lib['Tamborra_2014'] = 0.009*u.s
__sn_model_def_time_step_lib['Walk_2018'] = 0.005*u.s
__sn_model_def_time_step_lib['Walk_2019'] = 0.005*u.s
__sn_model_def_time_step_lib['Warren_2020'] = 0.015*u.s
__sn_model_def_time_step_lib['Zha_2021'] = 0.009*u.s #0.0009, but this seems waaay too small. we want ~200-500 bins


def sn_model_default_time_step(modelname: str) -> u.s:
    if modelname in __sn_model_def_time_step_lib.keys():
        return __sn_model_def_time_step_lib[modelname]
    else:
        print("Model name doesn't have a default time step, so using 1s as default")
        return 1.0*u.s

'''
Zha_2021
'''
fname = './SNEWPY_models/Zha_2021/s16.dat'
snewpy_models['Zha_2021'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Zha_2021/s17.dat',
        './SNEWPY_models/Zha_2021/s18.dat',
        './SNEWPY_models/Zha_2021/s19.dat',
        './SNEWPY_models/Zha_2021/s20.dat',
        './SNEWPY_models/Zha_2021/s21.dat',
        './SNEWPY_models/Zha_2021/s22.dat',
        './SNEWPY_models/Zha_2021/s23.dat',
        './SNEWPY_models/Zha_2021/s24.dat',
        './SNEWPY_models/Zha_2021/s25.dat',
        './SNEWPY_models/Zha_2021/s26.dat',
        './SNEWPY_models/Zha_2021/s30.dat',
        './SNEWPY_models/Zha_2021/s33.dat'
    ],
    'Zha_2021',
    Zha_2021
    )

'''
Warren_2020
'''
fname = './SNEWPY_models/Warren_2020/stir_a1.23/stir_multimessenger_a1.23_m10.0.h5'
snewpy_models['Warren_2020'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Warren_2020/stir_a1.23/stir_multimessenger_a1.23_m25.0.h5',
        './SNEWPY_models/Warren_2020/stir_a1.27/stir_multimessenger_a1.27_m10.0.h5',
        './SNEWPY_models/Warren_2020/stir_a1.27/stir_multimessenger_a1.27_m25.0.h5',
        './SNEWPY_models/Warren_2020/stir_a1.25/stir_multimessenger_a1.25_m55.h5'
    ],
    'Warren_2020',
    Warren_2020
    )

'''
Walk_2018
'''
fname = './SNEWPY_models/Walk_2018/s15.0c_3D_nonrot_dir1'
snewpy_models['Walk_2018'] = SNEWPYModel(
    [fname],
    'Walk_2018',
    Walk_2018
    )

'''
Tamborra_2014
'''
fname = './SNEWPY_models/Tamborra_2014/s20.0c_3D_dir1'
snewpy_models['Tamborra_2014'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Tamborra_2014/s27.0c_3D_dir1'
    ],
    'Tamborra_2014',
    Tamborra_2014
    )

'''
Sukhbold_2015
'''
fname = './SNEWPY_models/Sukhbold_2015/sukhbold-SFHo-z9.6.fits'
snewpy_models['Sukhbold_2015'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Sukhbold_2015/sukhbold-SFHo-s27.0.fits',
        './SNEWPY_models/Sukhbold_2015/sukhbold-LS220-s27.0.fits',
        './SNEWPY_models/Sukhbold_2015/sukhbold-LS220-z9.6.fits'
    ],
    'Sukhbold_2015',
    Sukhbold_2015
    )

'''
Kuroda_2020
'''
# TODO: note that there are several different models in this set
fname = './SNEWPY_models/Kuroda_2020/LnuR00B00.dat'
snewpy_models['Kuroda_2020'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Kuroda_2020/LnuR10B12.dat',
        './SNEWPY_models/Kuroda_2020/LnuR10B13.dat'
    ],
    'Kuroda_2020',
    Kuroda_2020
    )

'''
Fornax_2021
'''
fname = './SNEWPY_models/Fornax_2021/lum_spec_12M_r10000_dat.h5'
snewpy_models['Fornax_2021'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Fornax_2021/lum_spec_13M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_14M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_15M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_16M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_17M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_18M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_19M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_20M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_21M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_22M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_23M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_25M_r10000_dat.h5',
        './SNEWPY_models/Fornax_2021/lum_spec_26M_r10000_dat.h5'
    ],
    'Fornax_2021',
    Fornax_2021
    )

'''
Nakazato_2013
'''
fname = './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits'
snewpy_models['Nakazato_2013'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Nakazato_2013/nakazato-LS220-BH-z0.004-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev200ms-s13.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev300ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev200ms-s13.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev300ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-BH-z0.004-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev200ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev100ms-s13.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev200ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev300ms-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s13.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev200ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev100ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev200ms-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev300ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev300ms-s13.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev100ms-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev200ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-togashi-BH-z0.004-s30.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev300ms-s20.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev100ms-s50.0.fits',
        './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.02-t_rev300ms-s13.0.fits'
    ],
    'Nakazato_2013',
    Nakazato_2013
    )

'''
Bollig_2016
'''
fname = './SNEWPY_models/Bollig_2016/s11.2c'
snewpy_models['Bollig_2016'] = SNEWPYModel(
    [
        fname,
        './SNEWPY_models/Bollig_2016/s27.0c'
    ],
    'Bollig_2016',
    Bollig_2016
    )

'''
Walk_2019
'''
fname = './SNEWPY_models/Walk_2019/s40.0c_3DBH_dir1'
snewpy_models['Walk_2019'] = SNEWPYModel(
    [fname],
    'Walk_2019',
    Bollig_2016
    )


'''
Fornax_2019
'''
# fname = './SNEWPY_models/Fornax_2019/lum_spec_10M.h5'
# #TODO: fix issue with 'missing 2 required positional arguments'
# snewpy_models['Fornax_2019'] = SNEWPYModel(
#     fname,
#     'Fornax_2019',
#     Fornax_2019(fname)
#     )