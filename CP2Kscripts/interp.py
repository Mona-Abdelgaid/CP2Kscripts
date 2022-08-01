import click
from CP2Kscripts.utils import interp

@click.command(name = 'interp',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('initial_coordinates', type=str, nargs=1)
@click.argument('final_coordinates', type=str, nargs=1)
@click.argument('number_images', type=float, nargs=1)
@click.argument('output_name', type=str, nargs=1)

def main(initial_coordinates,final_coordinates,number_images,output_name):
    #The function requires:
        #  (1) initial corrdinates in xyz format i.e. 'initial_state.xyz' -- string 
        #  (2) final coordinates in xyz format i.e. 'final_state.xyz' -- string 
        #  (3) number of images to be generated between the inital and final coordinates i.e. 8 -- float 
        #  (4) output images name i.e. 'TS' -- string 

    interp(initial_coordinates,final_coordinates,number_images,output_name)