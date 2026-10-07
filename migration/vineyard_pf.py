import numpy as np
import yaml 
import matplotlib.pyplot as plt

# need to parse the phonopy qpoints.yaml files to extract frequencies at chosen q point/s
# yamls are convenient to parse at least!
# will need to disregard the soft negative modes found in some of these configurations...
# either due to not relaxing to low enough force tol, some inherent instability (sitting at higher min of the anharmonic PES), or finite size effects

# don't need this, dyn mat eigs should already be in THz
#conversion_factor_to_THz = 15.633302

data_min = yaml.safe_load(open("./min_left/qpoints.yaml"))
data_sad = yaml.safe_load(open("./saddle/qpoints.yaml"))
#dynmat_min = []
#dynmat_data_min = data_min["phonon"][0]["dynamical_matrix"]
#for row in dynmat_data_min:
#    vals = np.reshape(row, (-1, 2))
#    dynmat_min.append(vals[:, 0] + vals[:, 1] * 1j)
#dynmat_min = np.array(dynmat_min)

#(
#    eigvals_min,
#    eigvecs_min,
#) = np.linalg.eigh(dynmat_min)
#frequencies_min = np.sqrt(np.abs(eigvals_min.real))*np.sign(eigvals_min.real)

freq_data_min = data_min["phonon"][0]["band"]
freq_data_sad = data_sad["phonon"][0]["band"]
freq_min = np.zeros(len(freq_data_min))
freq_sad = np.zeros(len(freq_data_min))

#freq_prod_min = 1
#for j in range(len(freq_data_min)):
#    freq = freq_data_min[j]["frequency"]
#    if (j > 0) & (freq > 0):
#        freq_min[j] = freq
#    else:
#        freq_min[j] = 1
#        
#    freq_prod_min = freq_prod_min*freq_min[j]

def freq_prod(N: int, freq_data):
    prod = 1
    for j in range(N):
        freq = freq_data[j]["frequency"]
        if (j > 0) & (freq > 0):
            prod = prod*freq
        else:
            prod = prod
        
    return prod

omega_min = freq_prod(len(freq_data_min), freq_data_min)
omega_sad = freq_prod(len(freq_data_sad), freq_data_sad)
gamma0 = omega_min/omega_sad

#print(freq_min)
print(omega_min)
print(omega_sad)
print(gamma0)
