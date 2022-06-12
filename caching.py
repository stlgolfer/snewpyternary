#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:31:22 2022
Interal caching functions for snewpyternary
@author: phyics
"""

import os
import json
import numpy as np
import shutil

lock_file_template = {
    'fnames': []
    }

def load_lock_file():
    '''
    Loads the lock file and creates one if there isn't one available

    Returns
    -------
    dict
        The lock file as a dict

    '''
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
    '''
    Overwrite the lock file and store the information in JSON format

    Parameters
    ----------
    data : dict
        The lock file information to be locked

    Returns
    -------
    None.

    '''
    load_lock_file() # just makes sure the lock infrastructure is there
    with open('./.st_cache/cache.lock','w') as lock_file:
        lock_file.write(json.dumps(data))
        
def cache(fname, data):
    '''
    Caches data using a filename

    Parameters
    ----------
    fname : str
        The identifier to be used in the cache that coresponds to the real file
    data : numpy array
        Muse be a numpy array

    Returns
    -------
    None.

    '''
    # check if there is a lock file
    lock_file = load_lock_file()
    path = f'./.st_cache/{fname}'
    np.save(path,data)
    lock_file['fnames'].append(fname)
    # with open(path,'w') as file:
        
    # update lock file
    lock(lock_file)

def in_cache(name):
    '''
    Loads the lock file and tests whether a filename is in the lock file

    Parameters
    ----------
    name : str
        The fname to be checked

    Returns
    -------
    bool
        Whether the filename is in the lock file or not

    '''
    lock_file = load_lock_file()
    return name in lock_file['fnames']

def load_cache(fname):
    '''
    Loads a numpy array from the cache

    Parameters
    ----------
    fname : str
        The fname of the cached file

    Returns
    -------
    numpy array
        The numpy data from the cache

    '''
    lock_file = load_lock_file()
    if in_cache(fname):
        return np.load(f'./.st_cache/{fname}.npy',allow_pickle=True)
    else:
        # data was not available
        return None

def delete_cache():
    '''
    Completely removes the cache folder and all its contents

    Returns
    -------
    None.

    '''
    shutil.rmtree('./.st_cache')