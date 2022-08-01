import click
from CP2Kscripts.utils import gibbs_molecules

@click.command(name = 'gibbs_molecules',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('mol_file', type=str, nargs=1)


def main(mol_file):

    gibbs_molecules(mol_file)