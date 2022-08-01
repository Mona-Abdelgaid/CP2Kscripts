#!/usr/bin/env python

import argparse
import numpy as np
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
from sys import argv
import math
import ase.io
import ase
from ase import atoms
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
import ase.units
#=======================================================================
def interp(initial_coordinates,final_coordinates,number_images,output_name):
    #The function requires:
    #  (1) initial corrdinates in xyz format i.e. 'initial_state.xyz' -- string 
    #  (2) final coordinates in xyz format i.e. 'final_state.xyz' -- string 
    #  (3) number of images to be generated between the inital and final coordinates i.e. 8 -- float 
    #  (4) output images name i.e. 'TS' -- string 
    
    # Create ASE Atoms objects
    initial = ase.io.read(initial_coordinates, format="xyz")
    final = ase.io.read(final_coordinates, format="xyz")

    # Interpolate the position an end vectors
    interp_positions = np.linspace(initial.get_positions(), final.get_positions(), num=number_images)

    # Perform the interpolation
    images = []
    for position in interp_positions:
        image_atoms = initial.copy() # Copy to a new atoms object
        image_atoms.positions = position # Set the new atoms object's position
        images.append(image_atoms.copy()) # Add this new object to the list of images

    # Write the interpolated set of images to disk
    for number, image in enumerate(images):
        output_filename = f"{output_name}_{number}.xyz" # e.g. outputname_1, outputname_2, ..., etc.
        output_coordinates = ase.io.write(output_filename, image)
    return output_coordinates

#=======================================================================

def shift_coord(mol_file,title,out_file_name,vib_line,number_r_atoms,number_total_atoms,dispfac):
    # This function shifts the coordinates by a certain fraction in a certain direction, aiming to distrubt the system to remove a certain imaginary frequency 
    # 
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

    bohr2angs = 0.529177			#the conversion factor from bohr to angs

    with open(mol_file, 'r') as a: 		
        disp_x = []
        disp_y = []
        disp_z = []
        coord_element = []
        coord_x = []
        coord_y= []
        coord_z= []

        for line in a:
            if 'Vibration %d'%(vib_line) == line.strip() or  'vibration      %d'%(vib_line) == line.strip():
                for __ in range(number_r_atoms):
                    disp = next(a)			#reads the next lines in a after the line vibration %d
                    disp_data = disp.split()		#splits the displacement data into seperate items
                    disp_x.append(float(disp_data[0]))   #disp1, disp2, disp3 are the storage areas in which the displacement data from each column are put into
                    disp_y.append(float(disp_data[1]))
                    disp_z.append(float(disp_data[2]))
                    #print(disp_x, disp_y, disp_z)

                if number_r_atoms < number_total_atoms:
                    disp_x = disp_x + [0 for __ in range(number_total_atoms - number_r_atoms)]
                    disp_y = disp_y + [0 for __ in range(number_total_atoms - number_r_atoms)]
                    disp_z = disp_z + [0 for __ in range(number_total_atoms - number_r_atoms)]
        a.seek(0)

        for line in a:			
            if 'FR-COORD' in line:
                for __ in range(number_total_atoms):   	
                    coordinates = next(a)		#reads the next lines in a after the FR-COORD, coord in bohrs
                    coord_data = coordinates.split()	#splits the coordinates data into seperate items
                    coord_element.append(coord_data[0])		#coord1, coord2, coord3 are the storage areas in which the coordinates data from each column are put into
                    coord_x.append(float(coord_data[1]))
                    coord_y.append(float(coord_data[2]))
                    coord_z.append(float(coord_data[3]))
                    #print(coord_element, coord_x, coord_y, coord_z)

        disp_x = [i * dispfac for i in disp_x] #multiplies the displacements by a displacement factor 
        disp_y = [i * dispfac for i in disp_y]
        disp_z = [i * dispfac for i in disp_z]

        coord_x_new = [disp_x[i] + coord_x[i] for i in range(number_total_atoms)] #new coordinates in bohrs with the displacement factor added
        coord_y_new = [disp_y[i] + coord_y[i] for i in range(number_total_atoms)]
        coord_z_new = [disp_z[i] + coord_z[i] for i in range(number_total_atoms)]

        coord_x_new = ["{0:.10f}".format(i * bohr2angs) for i in coord_x_new] #new coordinates from bohrs to angstroms
        coord_y_new = ["{0:.10f}".format(i * bohr2angs) for i in coord_y_new]
        coord_z_new = ["{0:.10f}".format(i * bohr2angs) for i in coord_z_new]

    
    with open(out_file_name, 'w') as b:          #creates the xyz file
        b.write(str(number_total_atoms))
        b.write('\n')
        b.write(title)
        b.write('\n')
        for x in range(number_total_atoms):
            xyz = '{:<4}  {:>14}  {:>14}  {:>14}'.format(coord_element[x], coord_x_new[x], coord_y_new[x], coord_z_new[x])
            shifted_coordinate = b.write(xyz)
            shifted_coordinate = b.write('\n')
    
    return shifted_coordinate

#=======================================================================

def gibbs_molecules(mol_file):
    print('Molecules present in file:\n')
    molecules = ('co', 'co2', 'h2o', 'h2', 'carbonic_acid','bicarbonate', 'hydrodium', 'methylthiol', 'HSCH3', 'methane', 'hcl', 'ph3', 'HEthPh', 'propane','propylene')
    for count, molecules in enumerate(molecules):
        print(count,':', molecules)
    mol_type = int(input('Give the molecule number: '))
    
    #Data for molecules
    m_mol = [27.99491, 43.98983, 18.01056, 2.01565, 62.00039, 60.99257, 19.01839, 48.00337, 46.99555, 16.03130, 35.97668, 33.99724, 138.05032, 44.06260, 42.04695] #Mlecular_mass from gaussian units in g/mol
    m_mult = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,1]					#no units for multiplicity
    rot_cont = [2.67905, 0.54424, [36.75422, 20.81451, 13.28883], 87.72801, [ 0.55367, 0.53572, 0.27227], [0.57503, 0.53622, 0.27747], 
    [ 15.53711, 15.53614, 8.92453], [4.85477, 0.60872, 0.58359], [7.64508, 0.60703, 0.60702], [7.48285, 7.48280, 7.48279],14.82210, [6.22901, 6.22689, 5.59547],[0.20717, 0.02749, 0.02546], [1.40839, 0.40791, 0.36034],[2.24372, 0.44877, 0.39283]]			#from gaussian rotational constants in K
    num_rot = [ 1, 1, 3, 1, 3, 3, 3, 3, 3, 3, 1, 3,3,3,3]	#number of rotational constants
    sym_num = [ 1, 2, 2, 2, 2, 1, 3, 3, 3, 12, 1,3,1,2,6]      #symmetry number found online
    fu = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]  #fugacity, only listed number for water, rest were left as 0, currently, only water applies in calculations being in liquid phase
    num_vib = int(input('Number of vibrations:'))	#numer of vibrations
    na = 6.02214086E23				#avogadros constant mol^-1 
    m_kg =(m_mol[mol_type]/na) * (1.0/1000)		# mass,  grams

    #Relevant constants for calculations
    pi_num = 3.14159265359				# pi constant
    kb = 8.61733E-05				#boltzman constant, eV/K
    kb_d = 1.38064852E-23				#boltzman constant, (m^2 kg)/(s^2 * K)
    temp = 298					#temperature in kelvin
    pl_const = 4.135667662E-15			#plancks constant eV.s
    pl_const_d = 6.62607004E-34			#plancks constant, (m^2 kg) / s
    s_of_light = 2.998E10				#speed of light cm/s
    s_p = 101325					#standard pressure, Pa (kg/(m * s^2))
    pl_s = pl_const * s_of_light			#plancks constant * speed of light, eV.cm 

    #VIBRATIONAL CONTRIBUTIONS
    with open(mol_file, 'r') as a:
        vib = []
        for line in a:
            if '[FREQ]' in line:
                for x in range(num_vib):
                    vibs = next(a)                      #reads the next lines in a after the line vibration
                    vibs_data = vibs.split()
                    vib.append(float(vibs_data[0]))

    zpe_each = [0.5 * i * pl_s for i in vib]	#calculate the ZPE, eV
    zpe = sum(zpe_each)
    #heat capacity
    cv_each = [(i * pl_s) / (math.exp((i * pl_s)/(kb * temp))-1) for i in vib] 	#integrated heat capacity due to vibrational motion, eV
    cv = sum(cv_each)
    #entropy
    sv_each = [((i * pl_s) / ((kb * temp) * (math.exp((i * pl_s)/(kb * temp))-1)))-math.log(1-math.exp((-i * pl_s)/(kb * temp))) for i in vib]          #calculate sv in eV
    sv = kb * sum(sv_each)
    
    #TRANSLATIONAL CONTRIBUTIONS
    #heat capacity
    ct = (1.5) * kb * temp 		#integrated heat capacity due to translational motion
    #entropy
    qt = ((((2 * pi_num * m_kg * kb_d * temp) / (pl_const_d ** 2))) ** (1.5)) * ((kb_d * temp)/s_p) #translational partition function
    st = kb * (math.log(qt) + 5/2)      #entropy due to translational motion

    #ELECTRONIC CONTRIBUTIONS
    #heat capacity
    ce = 0 					#integrated heat capacity due to electronic motion
    #entropy
    qe = m_mult[mol_type]                        #electronic partition function, multiplicity
    se = kb * math.log(qe)              #entropy due to electronic motion also written se = kb * ln(2S+1)
    
    #ROTATIONAL CONTRIBUTIONS
    #heat capacity and entropy combined in if statement
    if num_rot[mol_type] == 1: 
        rota = rot_cont[mol_type]
        sym = sym_num[mol_type]
        qr = temp/(rota * sym)
        sr = kb * (math.log(qr) + 1)			#entropy due to rotational motion linear molecule
        cr = kb * temp					#integrated heat capacity due to rotational motion linear molecule

    elif num_rot[mol_type] == 3:
        rota = rot_cont[mol_type]
        sym = sym_num[mol_type]
        rot = 1
        for num in rot_cont[mol_type]:
            rot *= num
        qr = (math.sqrt(pi_num)/sym) * (((temp ** 3)/(rot)) ** 0.5)
        sr = kb * (math.log(qr) + (1.5))		#entropy due to rotational motion nonlinear molecule
        cr = (1.5) * kb * temp                         #integrated heat capacity due to rotational motion nonlinear molecule
        print ("rota", rota)
        print ("sym", sym)
        print ("sym", sym)
    else:
        print("check that rot_const is 1 or 3")

#For the free energy of a liquid substance
    if fu[mol_type] != 0:
        liq = math.log(fu[mol_type]/s_p)
    else: 
        liq = 0

    kbT = kb * temp 			#integrated boltzman to convert from cv to cp
    cp =  kbT + ct + cr + cv + ce		#heat capacity at constant pressure
    s_all = st + sr + se + sv - kb * liq    #summation of entropy

    t_s_all = temp * s_all			#TS contribution to free enery

    print ("kbT:", kbT)
    print ("ct:", ct)
    print ("cr:", cr)
    print ("cv:", cv)
    print ("ce:", ce)
    print ("st:", st)
    print ("sr:", sr)
    print ("sv:", sv)
    print ("se:", se)
    print ("The zpe, cp, and ts are:" , zpe, cp, t_s_all, "eV")	#print zpe, heat capacity, and entropy

#=======================================================================

def vib_list(out_file):
    # This function takes the .out file and extracts the vibrational modes and saves them in a .txt file -- 
    # The .txt file will be used in "thermo" script to get free energy corrections 
    with open(out_file, 'r') as f:
        vib = []
        for line in f:
            if 'VIB|Frequency (cm^-1)' in line:
                data = line.split()  #this splits all the lines that contains all the vibrational frequencies
                vib.append(float(data[2]))   #The following three lines stores a list of all freqs in vib
                vib.append(float(data[3]))
                vib.append(float(data[4]))

    FREQ = open('vib.txt', 'ab') #saving the frequencies in a text file
    np.savetxt(FREQ, vib, fmt = '%10.6f')
    FREQ.close()

#=======================================================================

def thermo(vib_list,temp):
    # The thermo function requires vibrations list in .txt format (all real frequencies, including img. freq. will result in an error )
    # The thermo function is for adsorbates -- uses Harmonic-oscillator-approximation

    # Read in vibrational frequencies (assumes cm^-1 units)
    vib = np.loadtxt(vib_list)

    # convert from cm-1 to eV
    vib_energies = vib * ase.units.invcm

    # initialize Harmonic-oscillator-approximated thermochemistry object
    thermo = HarmonicThermo(vib_energies)

    # Calculate U, ZPE, Cv, S, and ST.
    # Note that when E_pot = 0, U = ZPE + Cv
    u = thermo.get_internal_energy(temp, verbose=False)
    zpe = thermo.get_ZPE_correction()
    cv = u - zpe
    s = thermo.get_entropy(temp, verbose=False)
    ts = s * temp

    # calculate the G correction term
    g_correction = zpe + cv - ts

    # convert the results into a dictionary
    results = {'E_ZPE': zpe,
            'Cv_harm (0->T)': cv,
            'U': u,
            'S_harm': s,
            'T*S': ts,
            f'G_correction @ {temp} K': g_correction}

    # print the results in a nice format
    center = 50
    line = '=' * center
    print(line)
    print(f'Harmonic Thermochem Results at T = {temp} K'.center(center))
    print(line)

    # iteratively print each result in table format
    for name in results:
        lbl = 'eV/K' if name == 'S_harm' else 'eV'
        space = center - len(name) - 5 
        val = results[name]
        print(f'{name}{val:>{space}.9f} {lbl}')
    print(line)
#=======================================================================
