#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:44:31 2022
Going to figure out the dimensions of the scatter data that ternary plots uses
@author: phyics
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.models import Nakazato_2013
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
import os
import ternary
import math

scale=100
figure, tax = ternary.figure(scale=scale)
tax.boundary(linewidth=2.0)
tax.gridlines(color="blue", multiple=math.floor(scale/10))
tax.set_title("Test graph")
# pretty surer data is organized as (bottom,left,right)
### TODO: make sure that data_files[1] actually points to something that can get the header
# tax.bottom_axis_label(tables.get(data_files[1]).get("header").split(" ")[2])
# tax.left_axis_label(tables.get(data_files[1]).get("header").split(" ")[3])
# tax.right_axis_label(tables.get(data_files[1]).get("header").split(" ")[4])

plotting_data = [(20,80,0)] # top,right,left;
tax.scatter(points=plotting_data, color='green',label='yuh')
tax.ticks(axis='lbr', linewidth=1, multiple=scale/5)
tax.clear_matplotlib_ticks()

tax.show()