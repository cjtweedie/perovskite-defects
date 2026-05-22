import numpy as np
import os
from ase import Atoms
from ase import Atom
from ase.io.vasp import read_vasp
from ase.io.vasp import write_vasp

# TURN THIS INTO SOMETHING THAT CAN TAKE IN A PRISTINE POSCAR
# AND CHURN OUT A POSCAR WITH DEFECTS
# AFTER SPECIFYING THE Pb ATOM OF THE OCTAHEDRON OF CHOICE
# WILL GENERATE INTERSTITIALS ALONG ALL OCTAHEDRAL EDGES

# take in vasp POSCAR as atoms input object
poscar = read_vasp("POSCAR.vasp")

# find midpoint on line connecting two lattice sites, with given position vectors pos1, pos2
# can be in cartesian or direct (relative) coordinates
# (TO MAKE STUFF WORK LATER, MAKE SURE THIS IS CARTESIAN COORDS)
def interstitial_pos(pos1: np.array, pos2: np.array):
    r1 = 0.5*np.abs(pos1[0]+pos2[0])
    r2 = 0.5*np.abs(pos1[1]+pos2[1])
    r3 = 0.5*np.abs(pos1[2]+pos2[2])
    return np.array([r1, r2, r3]) 


def halide_pos(atoms: Atoms, idx_Pb: int):
    
    # need to find number of Pb atoms so can loop through this many * 2 halide indices
    for i in range(len(atoms)):
        if list(atoms.symbols)[i] == 'Pb':
            N_Pb = i
            break

    # can provide list of all halide atoms to calculate distances from each Pb to each halide 
    idx_halides = list(range(N_Pb*2, len(atoms)))
    
    # each octahedron contains central Pb atom and 6 surrounding halide atoms
    # for given Pb atom, find the 6 nearest-neighbour atoms 
    dists = atoms.get_distances(idx_Pb, idx_halides, mic=True)
        
    # sort halide dists by min dist but keep track of the indices too in halides_sorted
    dists_t = list(zip(dists, idx_halides))
    dists_t = sorted(dists_t)[0:6]
    halides_sorted = [j[1] for j in dists_t]
    halide_atoms = atoms[halides_sorted]
    
    # outputs the (Cartesian) positions of each halide atom around the given octahedron
    return np.array(halide_atoms.get_positions())

# choose which Pb atom to centre the defects around
# try to scan through a range of Pb atoms with different local chemical environments
# make this a function so can just interface through an i/o file
N = 39 # this is ion index in VESTA, which starts from 1, so need to -1 for python
halides = halide_pos(poscar, N-1)

# now generate all the positions of interstitial atoms,
# by placing them halfway between each n-n pair of halides in the octahedron
# find first 8 edges from 4 n-n atoms to each apical (max/min b vals) atom
# find last 4 from  2 n-n atoms to each 4c atom (max/min c vals) atom
# just do it manually LMAO

max_pos_a = np.argmax(halides[:,0])
min_pos_a = np.argmin(halides[:,0])    
max_pos_b = np.argmax(halides[:,1])
min_pos_b = np.argmin(halides[:,1])
max_pos_c = np.argmax(halides[:,2])
min_pos_c = np.argmin(halides[:,2])

int1 = interstitial_pos(halides[max_pos_b], halides[max_pos_c])
int2 = interstitial_pos(halides[max_pos_b], halides[max_pos_a])
int3 = interstitial_pos(halides[max_pos_b], halides[min_pos_c])
int4 = interstitial_pos(halides[max_pos_b], halides[min_pos_a])
int5 = interstitial_pos(halides[min_pos_b], halides[min_pos_a])
int6 = interstitial_pos(halides[min_pos_b], halides[max_pos_c])
int7 = interstitial_pos(halides[min_pos_b], halides[max_pos_a])
int8 = interstitial_pos(halides[min_pos_b], halides[min_pos_c])
int9 = interstitial_pos(halides[min_pos_c], halides[min_pos_a])
int10 = interstitial_pos(halides[min_pos_a], halides[max_pos_c])
int11 = interstitial_pos(halides[max_pos_c], halides[max_pos_a])
int12 = interstitial_pos(halides[max_pos_a], halides[min_pos_c])

# now just append these positions to end of the pristine atoms object, write POSCAR
# here assuming I interstitial
# FIX UP SO THIS IS MORE EFFICIENT/CLEANER -> FOR LOOPS?
# ALSO FIX SO THAT IT PUTS I INT IN THE RIGHT PLACE
# NEED TO REORDER I INDICES FOR I INT
poscar_Ii = poscar.copy()
poscar_Ii.append(Atom('I'))
poscar_Ii.positions[-1] = int1
write_vasp(f"I_int_Pb{N}_1.vasp", poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int2
#write_vasp("I_int_Pb%d_2.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int3
#write_vasp("I_int_Pb%d_3.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int4
#write_vasp("I_int_Pb%d_4.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int5
#write_vasp("I_int_Pb%d_5.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int6
#write_vasp("I_int_Pb%d_6.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int7
#write_vasp("I_int_Pb%d_7.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int8
#write_vasp("I_int_Pb%d_8.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int9
#write_vasp("I_int_Pb%d_9.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int10
#write_vasp("I_int_Pb%d_10.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int11
#write_vasp("I_int_Pb%d_11.vasp" % N, poscar_Ii, direct = True)
#poscar_Ii.positions[-1] = int12
#write_vasp("I_int_Pb%d_12.vasp" % N, poscar_Ii, direct = True)

# same thing but for Br interstitial
poscar_Bri = poscar.copy()
poscar_Bri.append(Atom('Br'))
poscar_Bri.positions[-1] = int1
write_vasp(f"Br_int_Pb{N}_1.vasp", poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int2
#write_vasp("Br_int_Pb%d_2.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int3
#write_vasp("Br_int_Pb%d_3.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int4
#write_vasp("Br_int_Pb%d_4.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int5
#write_vasp("Br_int_Pb%d_5.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int6
#write_vasp("Br_int_Pb%d_6.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int7
#write_vasp("Br_int_Pb%d_7.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int8
#write_vasp("Br_int_Pb%d_8.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int9
#write_vasp("Br_int_Pb%d_9.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int10
#write_vasp("Br_int_Pb%d_10.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int11
#write_vasp("Br_int_Pb%d_11.vasp" % N, poscar_Bri, direct = True)
#poscar_Bri.positions[-1] = int12
#write_vasp("Br_int_Pb%d_12.vasp" % N, poscar_Bri, direct = True)

# move all these files to output directory to clean folders up
dir_path = os.getcwd()
os.replace(f"{dir_path}/Br_int_Pb{N}_1.vasp", f"{dir_path}/../poscar_out/Br_int_Pb{N}_1.vasp")