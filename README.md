# CP2K scripts
## Installation

```console
root@host:~$ pip install CP2Kscripts
```

## Examples

- **Creating images for NEB**

The interp function requires: inital_state.xyz coordinates, final_state.xyz coordinates, number of replicas, output replicas name

    ```python
    CP2Kscripts.interp('IS.xyz','FS.xyz',5,'TS')
    ```

The output of this function are 5 replicas, namely, TS0.xyz, TS1.xyz, TS2.xyz, TS3.xyz, TS4.xyz

- **Shifting coordinates to remove extra imaginary frequency**

The shift_coord function requires: .mol file, text in the title of the xyz output file i.e. number of atoms in the system, name of the output file, order of the freq to be removed, number of atoms that will be shifted, total number of atoms in the system, and displacement fraction

    ```python
    CP2Kscripts.shift_coord('VIB-frequencies.mol-1.mol','139','new_shifted.xyz',2,139,139,1)
    ```

The output file of this function is new_shifted.xyz 

- **Extracting the vibational modes from CP2K output file**

The vib_list function requires .out file and generates vib.txt file which contains all the vibrational frequency modes

    ```python
    CP2Kscripts.vib_list('cp2k.out')
    ```

- **Free energy corrections for catalyst/or molecule adsorbed on catalyst (with no minor imaginary frequencies)**

The thermo function requires vib.txt file and temperature

    ```python
    CP2Kscripts.thermo('vib.txt',273.15)
    ```

The output is a summary table of all thermodynamic properties, i.e. entropy, enthalpy, zero-point energy, gibbs free energy corrections

- **Free energy corrections for catalyst/or molecule adsorbed on catalyst (with no minor imaginary frequencies)**

The gibbs_surface_all_R function requires .mol file and temperature

    ```python
    CP2Kscripts.gibbs_surface_all_R('VIB-frequencies.mol-1.mol',273.15)
    ```

The output of this function is the zero-point energy, heat capacity, and entropy in eV 

- **Free energy corrections for  catalyst/or molecule adsorbed on catalyst (with minor imaginary frequencies)**

The gibbs_surface function requires vib.txt file and temperature
The vib.txt file can be generated using the 'vib_list' function -- imaginary frequencies should be removed manually from the .txt file

    ```python
    CP2Kscripts.gibbs_surface('vib.txt',273.15)
    ```

The output of this function is the zero-point energy, heat capacity, and entropy in eV 

- **Free energy corrections for molecules**

The gibbs_molecules function requires vib.mol file for a certain molecule

    ```python
    CP2Kscripts.gibbs_molecules('vib.mol')
    ```

The output of this function is the zero-point energy, heat capacity, and entropy in eV 

- **Local electronic density of states at a particular site of a crystal**

The ldos function requires .pdos file of a particular atom in a given system

    ```python
    CP2Kscripts.ldos('PDOS-ALPHA_k1-1.pdos')
    ```
    