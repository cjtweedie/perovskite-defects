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

# should make this its own function file, if using in multiple scripts
# here though actually just getting indices...could generalise though
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
    
    # outputs the sorted list of halide indices only
    # don't need positions as we are just deleting halides at certain atom indices
    return halides_sorted


# choose which Pb atom to centre the defects around
# try to scan through a range of Pb atoms with different local chemical environments
# make this a function so can just interface through an i/o file
N = 57 # this is ion index in VESTA, which starts from 1, so need to -1 for python

# positions of all the halide vacancies is just the position vectors saved in halides
# NOPE here just indices
halides = halide_pos(poscar, N-1)
#print(halides[0])
#print(poscar.positions[halides[0]-1])
#print(poscar.positions[halides[0]])
#print(poscar.positions[halides[0]+1])

# don't know what species of halide X about to be removed, placeholder X for now
poscar_VX = poscar.copy()

# now for each halide position, need to remove that atom from the list to create vacancy
# need more systematic way so that no matter what atom index is, avoid global index mismatch
for i in range(len(halides)):    
    
    # need to identify if I/Br being removed, just get chem symbols from atoms object before deleting the atom
    # also figure out the Wyckoff position (if it was pristine cell) for nomenclature sake?
    V_sym = poscar_VX.symbols[halides[i]]
    #print(V_sym)
    del poscar_VX[halides[i]]
    
    # NEED TO REORDER INDICES FOR END POINTS SO THEY ARE ADJACENT TO START POINT INDICES
    # also figure out better naming scheme which automatically deals with degenerate halide positions from different Pb atoms
    write_vasp(f"{V_sym}_vac_Pb{N}_{halides[i]+1}.vasp", poscar_VX, direct = True)
    
    # move all these files to output directory as they are written to clean folders up
    os.replace(f"{V_sym}_vac_Pb{N}_{halides[i]+1}.vasp", f"{dir_path}/../poscar_out/{V_sym}_vac_Pb{N}_{halides[i]+1}.vasp")