import numpy as np
import yaml 
import matplotlib.pyplot as plt

# need to parse the phonopy qpoints.yaml files to extract frequencies at chosen q point/s
# yamls are convenient to parse at least!
# will need to disregard the soft negative modes found in some of these configurations...
# either due to not relaxing to low enough force tol, some inherent instability (sitting at higher min of the anharmonic PES), or finite size effects

# don't need this, dyn mat eigs should already be in THz
#conversion_factor_to_THz = 15.633302

data_min = yaml.safe_load(open("./min_left/qpoints_X.yaml"))
data_sad = yaml.safe_load(open("./saddle/qpoints_X.yaml"))

freq_data_min = data_min["phonon"][0]["band"]
freq_data_sad = data_sad["phonon"][0]["band"]

def freq_prod(N: int, freq_data):
    prod = 1
    for j in range(N):
        freq = freq_data[j]["frequency"]
        if freq > 0:
            freq = freq
        else:
            freq = 1
        
        prod = prod*freq
        
    return prod

omega_min = freq_prod(len(freq_data_min), freq_data_min)
omega_sad = freq_prod(len(freq_data_sad), freq_data_sad)
gamma0 = omega_min/omega_sad

print("omega_min =", omega_min)
print("omega_saddle =", omega_sad)
print("Vineyard prefactor =", gamma0)