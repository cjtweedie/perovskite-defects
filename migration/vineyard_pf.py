import numpy as np
import yaml 
import matplotlib.pyplot as plt

# need to parse the phonopy qpoints.yaml files to extract frequencies at chosen q point/s
# yamls are convenient to parse at least!
# will need to disregard the soft negative modes found in some of these configurations...
# either due to not relaxing to low enough force tol, some inherent instability (sitting at higher min of the anharmonic PES), or finite size effects

# don't need this, dyn mat eigs should already be in THz
#conversion_factor_to_THz = 15.633302

data_min = yaml.safe_load(open("./min/qpoints.yaml"))
data_sad = yaml.safe_load(open("./saddle/qpoints.yaml"))

freq_data_min = data_min["phonon"][0]["band"]
freq_data_sad = data_sad["phonon"][0]["band"]

# set some NONZERO threshold frequency below which phonons are excluded from the products
# ensures that acoustic modes at gamma are omitted even with numerical error
ftol = 0.05

def freq_scan(N_stop, freq_data):
    #scan through first few entries to remove the negative & acoustic modes
    for j in range(N_stop):
        freq = freq_data[j]["frequency"]
        if freq > ftol:
            break
        
    return j

# function to calculate the product of all non-negative phonon frequencies in the minimum and transition state configurations
def freq_prod(N: int, freq_data):
    N_start = freq_scan(5, freq_data)
    prod = 1
    for j in range(N_start, N):
        freq = freq_data[j]["frequency"]
        prod = prod * freq
        
    return prod

omega_min = freq_prod(len(freq_data_min), freq_data_min)
omega_sad = freq_prod(len(freq_data_sad), freq_data_sad)
pf = omega_min / omega_sad

f = open("vineyard_pf.txt", 'w')
print("frequency product (min) =", omega_min, "THz", file=f)
print("frequency product (saddle) =", omega_sad, "THz", file=f)
print("Vineyard prefactor =", pf, "THz", file=f)