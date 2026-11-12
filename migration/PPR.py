import numpy as np
import yaml 
import matplotlib.pyplot as plt

# READ FROM BANDS...NOT QPOINTS...NEED MASS DATA!
with open("./band.yaml", 'r') as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)
    
data_bands = data["phonon"][0]["band"]
data_atoms = data["points"]


# first, extract eigenvectors for given band index j, qpoint q (always Gamma here)
# this corresponds to phonon eigenmode e_jq(omega) 
def get_evec(phonon_bands, band_idx):
    N_atoms = len(phonon_bands[0]["eigenvector"])
    evec = [None]*N_atoms
    for k in range(N_atoms):
        evec[k] = [item[0] for item in(phonon_bands[band_idx]["eigenvector"][k])]
    return evec


def get_mass(atoms):
    N_atoms = len(atoms)
    mass_vec = [None]*N_atoms
    for k in range(N_atoms):
        mass_vec[k] = atoms[k]["mass"]
    return mass_vec


# calculate mod squared of atomic displacement vector within given normal mode e_jq
def mod_sq(evec):
    sum = 0
    for j in evec:
        sum = sum + j*j
    return sum


# functions for 1st & 2nd energy moments
# sum up the atom-projected kinetic energy terms for the e_jq eigenmode
# input mass vector m & displacement vector r, for all k 
# eigenvector contains r at every k as list of lists
# for the 1st energy moment, input m_k and |r_k|^2 as is, then square output
def M1(mass, disp):
    N_atoms = len(disp)
    sum = 0
    for k in range(N_atoms):
        sum = sum + mass[k]*mod_sq(disp[k])
    return sum*sum


# for the 2nd energy moment, use squares of the inputs
def M2(mass, disp):
    N_atoms = len(disp)
    sum = 0
    for k in range(N_atoms):
        sum = sum + (mass[k]*mod_sq(disp[k]))**2
    return sum  


e0 = get_evec(data_bands, 0)
m = get_mass(data_atoms)
print(e0[0])
print(m[0])
#print(mod_sq(e0[0]))

p0 = M1(m, e0)/(len(e0)*M2(m, e0))
print("M1 for band", 0, ":", M1(m, e0))
print("M2 for band", 0, ":", M2(m, e0))
print("PPR for band", 0, "at q = 0:", p0)