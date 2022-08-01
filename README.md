# CP2K scripts
## Installation

```console
root@host:~$ pip install CP2Kscripts
```

- **Creating images for NEB**
The interp function requires: inital_state.xyz coordinates, final_state.xyz coordinates, number of replicas, output replicas name

    ```python
    CP2Kscripts.interp('IS.xyz','FS.xyz',5,'TS')
    ```

- **Shifting coordinates to remove extra imaginary frequency**
The shift_coord function requires: vib.mol file, text in the title of the output file i.e. number of atoms in your system, name of the output file, order of the freq to be removed, number of atoms that will be shifted, total number of atoms in the system, and displacement fraction

    ```python
    CP2Kscripts.shift_coord('VIB-frequencies.mol-1.mol','139','new_shifted.xyz',2,139,139,1)
    ```

- **Extracting the vibational modes from CP2K output file**
The vib_list function requires .out file and generates vib.txt file which contains all the vibrational frequency modes

    ```python
    CP2Kscripts.vib_list('cp2k.out')
    ```
- **Free energy corrictions for catalyst/or molecule adsorbed on catalyst (with no minor imaginary frequencies)**
The thermo functionrequires vib.txt file and temperature

    ```python
    CP2Kscripts.thermo('vib.txt',273.15)
    ```
- **Free energy corrictions for molecules**
The gibbs_molecules function requires vib.mol file for a certain molecule

    ```python
    CP2Kscripts.gibbs_molecules('vib.mol')
    ```
- **Local electronic density of states at a particular site of a crystal**
The ldos function requires .pdos file of a particular atom in a given system

    ```python
    CP2Kscripts.ldos('PDOS-ALPHA_k1-1.pdos')
    ```