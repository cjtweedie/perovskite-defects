import numpy as np
import os
from ase import Atoms
from ase import Atom
from ase.io.vasp import read_vasp
from ase.io.vasp import write_vasp

# directory path variable for os stuff later
dir_path = os.getcwd()

# take in vasp POSCAR as atoms input object
# might need to change ASE read vasp behaviour so that it can read files from arbitrary directories
#poscar = read_vasp(f"{dir_path}/../poscar_in/POSCAR.vasp")
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
N = 57 # this is ion index in VESTA, which starts from 1, so need to -1 for python
halides = halide_pos(poscar, N-1)

# now generate all the positions of interstitial atoms,
# by placing them halfway between each n-n pair of halides in the octahedron
# just do it manually LMAO
max_pos_a = np.argmax(halides[:,0])
min_pos_a = np.argmin(halides[:,0])    
max_pos_b = np.argmax(halides[:,1])
min_pos_b = np.argmin(halides[:,1])
max_pos_c = np.argmax(halides[:,2])
min_pos_c = np.argmin(halides[:,2])

# assign these position vectors an index 0->11 in dictionary so they can be referenced in loops lateer
int_positions = {
    0: interstitial_pos(halides[max_pos_b], halides[max_pos_c]),
    1: interstitial_pos(halides[max_pos_b], halides[max_pos_a]),
    2: interstitial_pos(halides[max_pos_b], halides[min_pos_c]),
    3: interstitial_pos(halides[max_pos_b], halides[min_pos_a]),
    4: interstitial_pos(halides[min_pos_b], halides[min_pos_a]),
    5: interstitial_pos(halides[min_pos_b], halides[max_pos_c]),
    6: interstitial_pos(halides[min_pos_b], halides[max_pos_a]),
    7: interstitial_pos(halides[min_pos_b], halides[min_pos_c]),
    8: interstitial_pos(halides[min_pos_c], halides[min_pos_a]),
    9: interstitial_pos(halides[min_pos_a], halides[max_pos_c]),
    10: interstitial_pos(halides[max_pos_c], halides[max_pos_a]),
    11: interstitial_pos(halides[max_pos_a], halides[min_pos_c])
    }

# now just append these positions to end of the pristine atoms object, write POSCAR
# reference the corresponding position vector for each edge index 0->11, easily write unique file names in each loop with f strings
# do for both I and Br as the interstitial species
poscar_Ii = poscar.copy()
poscar_Bri = poscar.copy()
poscar_Ii.append(Atom('I'))
poscar_Bri.append(Atom('Br'))
for i in range(12):
    poscar_Ii.positions[-1] = int_positions[i]
    poscar_Bri.positions[-1] = int_positions[i]
    
    # ALSO FIX SO THAT IT PUTS I INT IN THE RIGHT PLACE
    # NEED TO REORDER I INDICES FOR I INT
    write_vasp(f"I_int_Pb{N}_{i+1}.vasp", poscar_Ii, direct = True)
    write_vasp(f"Br_int_Pb{N}_{i+1}.vasp", poscar_Bri, direct = True)
    
    # move all these files to output directory as they are written to clean folders up
    os.replace(f"{dir_path}/I_int_Pb{N}_{i+1}.vasp", f"{dir_path}/../poscar_out/I_int_Pb{N}_{i+1}.vasp")
    os.replace(f"{dir_path}/Br_int_Pb{N}_{i+1}.vasp", f"{dir_path}/../poscar_out/Br_int_Pb{N}_{i+1}.vasp")