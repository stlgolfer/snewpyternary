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

# snewpy-snowglobes stuff
import snewpy.snowglobes as snowglobes

modelFilePathBase = "./SNEWPY_models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
snowglobes_out_name = "snowglobes-output"
model

# we have the model loaded, so let's try to generate the fluence files using SNEWPY
snowglobes.generate_fluence(model_path=modelFilePath,model_type="Nakazato_2013",
                            transformation_type="NoTransformation",
                            d=10,
                            output_filename=snowglobes_out_name
                            )
snowglobes.collate(tarball_path=modelFilePathBase+snowglobes_out_name,
                   smearing=False,skip_plots=False,SNOwGLoBESdir='$SNOWGLOBES')