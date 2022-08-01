import click
from CP2Kscripts.utils import vib_list

@click.command(name = 'vib_list',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('out_file', type=str, nargs=1)


def main(out_file):

    vib_list(out_file)
