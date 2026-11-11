import numpy as np
import yaml 
import matplotlib.pyplot as plt

# READ FROM BANDS...NOT QPOINTS...NEED MASS DATA!
data = yaml.safe_load(open("./band.yaml"))

data_bands = data["phonon"][0]["band"]
#print([item[0] for item in(data_bands[0]["eigenvector"][0])])
#print(len(data_bands))

# get the masses from each atom to normalise the eigenmode displacement amplitudes

# first, extract eigenvectors for given band index j, qpoint q (always Gamma here)
# this corresponds to phonon eigenmode e_jq(omega) 
def get_evec(phonon_bands, N_atoms, band_idx):
    evec = [None]*N_atoms
    for j in range(N_atoms):
        evec[j] = [item[0] for item in(phonon_bands[band_idx]["eigenvector"][j])]
    
    return evec

# calculate mod squared of atomic displacement vector within given normal mode e_jq
def mod_sq(evec):
    sum = 0
    for j in evec:
        sum = sum + j*j
    return sum

e0 = get_evec(data_bands, 159, 0)
print(e0)
print(mod_sq(e0[0]))