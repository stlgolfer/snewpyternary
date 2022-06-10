#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:31:22 2022
Used for testing/creating caching features
@author: phyics
"""

import os
import json
import numpy as np

lock_file_template = {
    'fnames': []
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
        
def cache(fname, data):
    # check if there is a lock file
    lock_file = load_lock_file()
    path = f'./.st_cache/{fname}'
    np.save(path,data)
    lock_file['fnames'].append(fname)
    # with open(path,'w') as file:
        
    # update lock file
    lock(lock_file)

def in_cache(name):
    lock_file = load_lock_file()
    return name in lock_file['fnames']

def load_cache(fname):
    lock_file = load_lock_file()
    if in_cache(fname):
        return np.load(f'./.st_cache/{fname}.npy',allow_pickle=True)
    else:
        # data was not available
        return None