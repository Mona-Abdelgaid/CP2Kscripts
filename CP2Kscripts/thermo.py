import click
from CP2Kscripts.utils import thermo

@click.command(name = 'thermo',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('vib_list', type=str, nargs=1)
@click.argument('temp', type=float, nargs=1)

def main(vib_list,temp):
    # The thermo function requires vibrations list in .txt format (all real frequencies, including img. freq. will result in an error )
    # The thermo function is for adsorbates -- uses Harmonic-oscillator-approximation
    
    thermo(vib_list,temp)