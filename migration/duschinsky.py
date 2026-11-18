import numpy as np
import yaml 
from yaml import CLoader as Loader
import matplotlib.pyplot as plt
import matplotlib.colorizer as clr

with open("./min/qpoints.yaml", 'r') as f:
    data_min = yaml.load(f, Loader=Loader)
with open("./saddle/qpoints.yaml", 'r') as f:
    data_sad = yaml.load(f, Loader=Loader)

bands_min = data_min["phonon"][0]["band"]
bands_sad = data_sad["phonon"][0]["band"]


# first, extract eigenvectors for given band index j, qpoint q (always Gamma here)
# this corresponds to phonon eigenmode e_jq(omega) 
# need to extract in one array for ease of manipulation later though
def get_evec(phonon_bands, band_idx):
    N_atoms = len(phonon_bands[0]["eigenvector"])
    evec = [None]*N_atoms*3
    count = 1
    for k in range(N_atoms):
        count = count - 1
        for j in range(3):
            i = k + count
            #ek = [item[0] for item in(phonon_bands[band_idx]["eigenvector"][k])]
            #evec[i] = ek[j]
            evec[i] = [item[0] for item in(phonon_bands[band_idx]["eigenvector"][k])][j]
            count = count + 1
    return evec

# dot product function; vectors must be same length
def dot_prod(v, u):
    N = len(v)
    sum = 0
    for j in range(N):
        sum = sum + v[j]*u[j]
    return sum

# loop over all min/saddle eigenvectors
# setting up Duschinsky matrix as flattened array of length Nj^2
# fill entries in single loop
# for colormap, need matrix...lol
Nj = len(bands_min)
J_mat = [None]*Nj*Nj
i = 0
for k in range(Nj*Nj):
    j = k%Nj
    if j == 0 and k != 0:
        i = i+1
    J_mat[k] = abs(dot_prod(get_evec(bands_min, i), get_evec(bands_sad, j)))
    
J_mat = np.transpose(np.array(J_mat).reshape(Nj, Nj))

# range of eigenmodes to plot
N = 11
X = np.arange(N)
Y = np.flip(X)
fig1, ax1 = plt.subplots(1, 1)
ax1.pcolormesh(X, X, J_mat[0:N, 0:N], cmap='viridis', vmin=np.min(J_mat), vmax=np.max(J_mat), shading='auto')
cr = clr.Colorizer(cmap='viridis')
cbar = fig1.colorbar(clr.ColorizingArtist(cr), ax=ax1)
ax1.set_xlim(0-0.5, N-1+0.5)
ax1.set_ylim(N-1+0.5, 0-0.5)
ax1.set_xticks(X)
ax1.set_yticks(Y)
ax1.set_xlabel("Minimum eigenmode", size="medium")
ax1.set_ylabel("Saddle eigenmode", size="medium")
ax1.set_title("Duschinsky matrix")
#cbar.set_label("Overlap", labelpad=20, size="large")
fig1.tight_layout()
fig1.savefig("./duschinsky_matrix.png")