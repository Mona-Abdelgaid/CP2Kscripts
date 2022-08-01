import click
from CP2Kscripts.utils import shift_coord

@click.command(name = 'shift_coord',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('mol_file', type=str, nargs=1)
@click.argument('title', type=str, nargs=1)
@click.argument('out_file_name', type=str, nargs=1)
@click.argument('vib_line', type=str, nargs=1)
@click.argument('number_r_atoms', type=float, nargs=1)
@click.argument('number_total_atoms', type=float, nargs=1)
@click.argument('dispfac', type=float, nargs=1)

def main(mol_file,title,out_file_name,vib_line,number_r_atoms,number_total_atoms,dispfac):
    # This function shifts the coordinates by a certain fraction in a certain direction, aiming to distrubt the system to remove a certain imaginary frequency 
    
    # mol_file is vibration mol file name i.e. 'VIB-frequencies.mol-1.mol'  -- string 
    # title is the title line in the xyz file which is the number of atoms in the system i.e. '139' -- string 
    # out_file_name is the name of the xyz file which will hold the new coordinates i.e. 'shifted_coordinates.xyz' -- string 
    # vib_line is the "Vibration number" for displacement -- vibration line that the displacements will be taken from 
        # i.e. 2 (second vib number in case of TS because we want to keep the first vib number but want to remove the extra imaginary number)
    # number_r_atoms is the total number of atoms that will be displaced, all or selected unfrozen -- the number of atoms relevant to displacement i.e. 139 -- float
        #  (I want to displace all atoms i
        # n my system, 
        #  if the surface is fixed then that means only adsorbates will be shifted since the surface atoms will have 0 0 0 in the xyz direction, i.e. will not be shifted)
    # number_total_atoms is the total number of atoms in system -- the number of atoms overall i.e. 139 -- float 
    # dispfac is the displacement factor i.e. 1 -- float 

    shift_coord(mol_file,title,out_file_name,vib_line,number_r_atoms,number_total_atoms,dispfac)