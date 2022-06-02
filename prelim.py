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

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
snowglobes_out_name = "snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
tball_suffix = '.tar.bz2'
#print(os.environ['SNOWGLOBES'])
smearing = False
model

# we have the model loaded, so let's try to generate the fluence files using SNEWPY
snowglobes.generate_fluence(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=10,
                            output_filename=snowglobes_out_name
                            )

# also need to run the simulation so that we get the correct output
snowglobes.simulate(SNOwGLoBESdir=snowglobes_dir,
                    tarball_path=modelFilePathBase+snowglobes_out_name+tball_suffix,
                    detector_effects=smearing,
                    detector_input="all"
                    )

snowglobes.collate(tarball_path=modelFilePathBase+snowglobes_out_name+tball_suffix,
                   smearing=False,skip_plots=False,SNOwGLoBESdir=snowglobes_dir)