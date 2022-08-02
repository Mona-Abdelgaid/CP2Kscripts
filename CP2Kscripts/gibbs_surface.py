import click
from CP2Kscripts.utils import gibbs_surface

@click.command(name = 'gibbs_surface_all_R',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('vib_file', type=str, nargs=1)
@click.argument('temp', type=float, nargs=1)

def main(vib_file,temp):

    gibbs_surface(vib_file,temp)