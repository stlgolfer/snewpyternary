#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:31:22 2022
Used for testing/creating caching features
@author: phyics
"""

import os
import json

lock_file_template = {
    'available': {
        'detector_cache':False,
        'flux_cache':False
        }
    }

def load_lock_file():
    if (os.path.isdir('.st_cache')):
        # cache exists, so load it
        with open('./.st_cache/cache.lock','r') as lock_file:
            return json.load(lock_file)
    else:
        # create cache dir
        os.makedirs('.st_cache')
        # create a new lock file
        with open('./.st_cache/cache.lock','w') as lock_file:
            lock_file.write(json.dumps(lock_file_template))
            return lock_file_template

def lock(data):
    load_lock_file() # just makes sure the lock infrastructure is there
    with open('./.st_cache/cache.lock','w') as lock_file:
        lock_file.write(json.dumps(data))
        
def cache(fname, plot_data,raw_data,l_data):
    # check if there is a lock file
    lock_file = load_lock_file()
    # overwrite current data and update lock file
    organized = {
        'plot_data': plot_data,
        'raw_data': raw_data,
        'l_data': l_data
        }
    with open(f'./.st_cache/{fname}.json','w') as detector_file:
        detector_file.write(json.dumps(organized))
    # update lock file
    lock_file['available'][f'{fname}_cache'] = True
    lock(lock_file)

def load_cache(fname):
    lock_file = load_lock_file()
    if lock_file['available'][f'{fname}_cache']==True:
        with open(f'./.st_cache/{fname}.json','r') as detector_file:
            return json.load(detector_file)
    else:
        # data was not available
        return None       
    
cache('detector',[0,4],[1,2],[3,4])
detector_cache = load_cache('detector')