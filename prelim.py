# goal of this preliminary chart will be to generate fluences given a model
# then create a ternary model of fluences and then also the event rates
# first, let's import the snewpy stuff
# let's start with the Nakazato model
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

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

# simulation details
modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
ntbins = 20 # number of time bins
deltat = (1*u.s) # time bin size, details found in SNEWPY article
d = 10 # in pc, distance to SN
detector = "halo1" # refrain from using "all"

snowglobes_out_name = "snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
tball_suffix = 'kpc.tar.bz2'
#print(os.environ['SNOWGLOBES'])
smearing = True
model

'''
# we have the model loaded, so let's try to generate the fluence files using SNEWPY
snowglobes.generate_fluence(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=10,
                            output_filename=snowglobes_out_name
                            )
'''

tball_complete = snowglobes.generate_time_series(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=d,
                            output_filename=snowglobes_out_name,
                            ntbins=ntbins,
                            deltat=deltat
                            )

# also need to run the simulation so that we get the correct output
snowglobes.simulate(SNOwGLoBESdir=snowglobes_dir,
                    tarball_path=tball_complete,
                    detector_effects=smearing,
                    detector_input=detector
                    )

'''
this method is too collated--as a first try we just need the event calculations
'''
tables = snowglobes.collate(tarball_path=modelFilePathBase+snowglobes_out_name+tball_suffix,
                   smearing=smearing,skip_plots=True,SNOwGLoBESdir=snowglobes_dir)

data_files = list(tables.keys())
plotting_data = []
scale = 100

'''
so what we have to do is go through each time bin, sum each particle's event
count, then append to plotting_data. Then, we will have a list of time-evolved
count
'''

for file in data_files: # which effectively goes through each time bin
    # we want to only have the ones that end in 'unweighted' for now
    if ("unweighted" in file):
        # now fetch the time bin's energy chart
        e_bins = tables.get(file)
        
        # first we need to get the points by iterating through the points
        header = list(e_bins.get("header").split(" ")) # header in dictionary
        #tables.get(file).data[]
        data = e_bins.get("data") # 2D array for data
        # now go through each available particles
        a=np.sum(data[1])
        b=np.sum(data[2])
        c=np.sum(data[3])
        total = a+b+c
        plotting_data.append((100*a/total,100*b/total,100*c/total))
        '''
        if (max(a,b,c)>scale):
            scale=math.ceil(max(a,b,c))
        
        for x in range(len(data[0])+1,3): # we can exclude the energy bins
            proto_tuple.append((np.sum(data[x])))
        plotting_data.append(tuple(proto_tuple)) # converts that into a tuple
        '''

# then we can generate a ternary plot for this
figure, tax = ternary.figure(scale=scale)
tax.boundary(linewidth=2.0)
tax.gridlines(color="blue", multiple=math.floor(scale/5))
tax.set_title("Nakazato_2013")
# data is organized in top, right, left

### TODO: make sure that data_files[1] actually points to something that can get the header
tax.bottom_axis_label(tables.get(data_files[1]).get("header").split(" ")[1])
tax.left_axis_label(tables.get(data_files[1]).get("header").split(" ")[2])
tax.right_axis_label(tables.get(data_files[1]).get("header").split(" ")[3])

tax.scatter(points=plotting_data, color='green',label='yuh')
tax.ticks(axis='lbr', linewidth=1, multiple=math.floor(scale/5))
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off') # disables regular matlab plot axes

tax.show()