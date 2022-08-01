import click
from CP2Kscripts.utils import ldos

@click.command(name = 'ldos',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file', type=str, nargs=1)


def main(file):

    ldos(file)
