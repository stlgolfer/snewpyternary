#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 23:43:12 2022

@author: physics
"""

from snewpy.models import Bollig_2016, Fornax_2019, Fornax_2021, Kuroda_2020, Nakazato_2013, Sukhbold_2015, Tamborra_2014, Walk_2018, Walk_2019, Warren_2020, Zha_2021

class SNEWPYModel:
    """
    SNEWPY models have several available configs, so it's convenient to just
    bundle this information togeter
    """
    
    def __init__(self,file_path,model_type,model):
        self.file_path = file_path
        self.model_type = model_type
        self.model = model

snewpy_models = {}

#TODO: must also include a time step that makes sense. let's use a seperate dictionary for this to keep things organized

'''
Zha_2021
'''
fname = './SNEWPY_models/Zha_2021/s17.dat'
snewpy_models['Zha_2021'] = SNEWPYModel(
    fname,
    'Zha_2021',
    Zha_2021(fname)
    )

'''
Warren_2020
'''
fname = './SNEWPY_models/Warren_2020/stir_a1.23/stir_multimessenger_a1.23_m10.0.h5'
snewpy_models['Warren_2020'] = SNEWPYModel(
    fname,
    'Warren_2020',
    Warren_2020(fname)
    )

'''
Walk_2018
'''
fname = './SNEWPY_models/Walk_2018/s15.0c_3D_nonrot_dir1'
snewpy_models['Walk_2018'] = SNEWPYModel(
    fname,
    'Walk_2018',
    Walk_2018(fname)
    )

'''
Tamborra_2014
'''
fname = './SNEWPY_models/Tamborra_2014/s20.0c_3D_dir1'
snewpy_models['Tamborra_2014'] = SNEWPYModel(
    fname,
    'Tamborra_2014',
    Tamborra_2014(fname)
    )

'''
Sukhbold_2015
'''
fname = './SNEWPY_models/Sukhbold_2015/sukhbold-SFHo-z9.6.fits'
snewpy_models['Sukhbold_2015'] = SNEWPYModel(
    fname,
    'Sukhbold_2015',
    Sukhbold_2015(fname)
    )

'''
Kuroda_2020
'''
# TODO: note that there are several different models in this set
fname = './SNEWPY_models/Kuroda_2020/LnuR00B00.dat'
snewpy_models['Kuroda_2020'] = SNEWPYModel(
    fname,
    'Kuroda_2020',
    Kuroda_2020(fname)
    )

'''
Fornax_2021
'''
fname = './SNEWPY_models/Fornax_2021/lum_spec_12M_r10000_dat.h5'
snewpy_models['Fornax_2021'] = SNEWPYModel(
    fname,
    'Fornax_2021',
    Fornax_2021(fname)
    )

'''
Nakazato_2013
'''
snewpy_models['Nakazato_2013'] = SNEWPYModel(
    './SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits',
    'Nakazato_2013',
    Nakazato_2013('./SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits')
    )

'''
Bollig_2016
'''
fname = './SNEWPY_models/Bollig_2016/s11.2c'
snewpy_models['Bollig_2016'] = SNEWPYModel(
    fname,
    'Bollig_2016',
    Bollig_2016(fname)
    )

'''
Walk_2019
'''
fname = './SNEWPY_models/Walk_2019/s40.0c_3DBH_dir1'
snewpy_models['Walk_2019'] = SNEWPYModel(
    fname,
    'Walk_2019',
    Bollig_2016(fname)
    )


'''
Fornax_2019
'''
fname = './SNEWPY_models/Fornax_2019/lum_spec_10M.h5'
snewpy_models['Fornax_2019'] = SNEWPYModel(
    fname,
    'Fornax_2019',
    Fornax_2019(fname)
    )