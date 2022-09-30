
from snewpy.flavor_transformation import *
import data_handlers as handlers
import snewpyternary as t
import os
import ternary
import snowglobes_wrapper
import click
from model_wrappers import snewpy_models, sn_model_default_time_step

@click.command()
@click.argument('models',required=True,type=str,nargs=-1)
@click.option('--setno', required=False, default=[0], multiple=True, help="Set numbers")
@click.option('--all', required=False, default=False, help="For all set numbers. Will override setno param")
def start(models, setno, all):
    for model in models:
        print(model)
        print("==================================================\n")
        to_iterate = setno if all==False else range(len(snewpy_models[model].file_paths))
        for s in to_iterate:

            print(f'{model}:{s} | {snewpy_models[model].file_paths[s]}')
            print(str(snewpy_models[model].model(snewpy_models[model].file_paths[s])))
            print("------------------------------------------------\n")

if __name__ == '__main__': # for multiprocessing
    start()