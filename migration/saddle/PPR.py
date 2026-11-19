import numpy as np
import yaml 
from yaml import CLoader as Loader, CDumper as Dumper
import matplotlib.pyplot as plt

# DON'T NEED MASSES ANYMORE, read from qpoints w eigenvectors on
with open("./qpoints.yaml", 'r') as f:
#    data = yaml.load(f, Loader=yaml.SafeLoader)
    data = yaml.load(f, Loader=Loader)
data_bands = data["phonon"][0]["band"]

# still need mass for some things...but read from phonopy_disp
with open("./phonopy_disp.yaml", 'r') as f:
    data = yaml.load(f, Loader=Loader)
data_atoms = data["primitive_cell"]["points"]

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
def mod_sq(vector):
    sum = 0
    for j in vector:
        sum = sum + j*j
    return sum


# functions for 1st & 2nd energy moments
# phonopy outputs eigenvectors already normalised; so masses not necessary
# input atom-projected components of the eigenvector = atomic displacement vector, for all k for e_jq eigenmode
# for the 1st moment, square output
def M1(disp):
    N_atoms = len(disp)
    sum = 0
    for k in range(N_atoms):
        #sum = sum + mass[k]*mod_sq(disp[k])
        sum = sum + mod_sq(disp[k])
    return sum*sum


# for the 2nd energy moment, use squares of the inputs
def M2(disp):
    N_atoms = len(disp)
    sum = 0
    for k in range(N_atoms):
        #sum = sum + (mass[k]*mod_sq(disp[k]))**2
        sum = sum + (mod_sq(disp[k]))*(mod_sq(disp[k]))
    return sum  

## calculate PPR ##

# very first mode is the unstable negative freq, highly localised on migrating ion
# wants to move from local max of PES (TS) to local min 
e0 = get_evec(data_bands, 0)

# acoustic modes [0,2] aren't relevant; at q=0 they should all be collective translations
# so PPR ~ 1
e1 = get_evec(data_bands, 1)

# first optical mode, probably more localised? so PPR < 1
e4 = get_evec(data_bands, 4)

p0 = M1(e0)/(len(e0)*M2(e0))
p1 = M1(e1)/(len(e1)*M2(e1))
p4 = M1(e4)/(len(e4)*M2(e4))

print("(sum of |n|^2)^2 (unstable) = ", M1(e0))
print("sum of (|n|^2)^2 (unstable) = ", M2(e0))

print("(sum of |n|^2)^2 (acoustic) = ", M1(e1))
print("sum of (|n|^2)^2 (acoustic) = ", M2(e1))

print("PPR for band 1 (unstable)", "at q = 0:", p0)
print("PPR for band 2 (acoustic)", "at q = 0:", p1)
print("PPR for band 5 (optical)", "at q = 0:", p4)

# check normalisation of eigenmodes; need to loop over all elements from all atomic contributions
e0_full = [None]*len(e0)*3
count = 1
for k in range(len(e0)):
    count = count -1
    for j in range(3):
        i = k + count
        e0_full[i] = e0[k][j]
        count = count + 1

print("|e0|^2 = ", mod_sq(e0_full))

# loop over all phonon modes j, plot the PPRs?
# finally make an Natoms variable, won't be too slow here because not referencing in different scope
Nj = len(data_bands)
Nk = len(data_bands[0]["eigenvector"])
freqs = [None]*Nj
modes = [None]*Nj
pj = [None]*Nj
for j in range(Nj):
    ej = get_evec(data_bands, j)
    pj[j] = M1(ej)/(Nk*M2(ej)) 
    freqs[j] = data_bands[j]["frequency"] 
    modes[j] = j+1
    
    
# also plot the magnitude of atom projections, for each mode
# each atomic displacement component r_k
# plot the projection into each of 159 atoms for all 477 modes??
# or show proj into atomic species Cs, Pb, I, Br
# since I only care about behaviour of individual modes, can look at all atoms
# also plot the magnitude of atom projections, for each mode
# each atomic displacement component r_k
# plot the projection into each of 159 atoms for all 477 modes??
# or show proj into atomic species Cs, Pb, I, Br
# since I only care about behaviour of individual modes, can look at all atoms
# need to normalise the displacement vectors by atomic mass
proj_e0 = [None]*Nk
proj_e1 = [None]*Nk
proj_e4 = [None]*Nk
atoms = [None]*Nk
m = get_mass(data_atoms)
f0 = open("proj_e0", 'w')
f1 = open("proj_e1", 'w')
f4 = open("proj_e4", 'w')
print("Atom no.", "Eigenmode projection weight", file=f0)
print("Atom no.", "Eigenmode projection weight", file=f1)
print("Atom no.", "Eigenmode projection weight", file=f4)
for k in range(Nk):
    proj_e0[k] = mod_sq(e0[k])/m[k]
    proj_e1[k] = mod_sq(e1[k])/m[k]
    proj_e4[k] = mod_sq(e4[k])/m[k]
    atoms[k] = k+1
    print(atoms[k], proj_e0[k], file=f0)
    print(atoms[k], proj_e1[k], file=f1)
    print(atoms[k], proj_e4[k], file=f4)
    
# mode range for PPR plot
N = 10
fig2, ax2 = plt.subplots(1, 1)
ax2.scatter(modes[0:N], pj[0:N], linewidth=1.0) 
ax2.set_xlabel("Phonon mode index j", size="medium")
ax2.set_ylabel("Phonon participation ratio", size="medium")
ax2.set_xticks(np.arange(N))
ax2.set_title("VI(0) 8d-I4c, min")
fig2.tight_layout()
fig2.savefig("./PPR_low10.png")

fig3, ax3 = plt.subplots(1, 1)
ax3.scatter(atoms, proj_e0, linewidth=1.0) 
ax3.set_xlabel("Atom index", size="medium")
ax3.set_ylabel("Atomic projection", size="medium")
ax3.set_title("VI(0) 8d-I4c, saddle, e0 (unstable mode)")
fig3.tight_layout()
fig3.savefig("./atom_proj_e0.png")

fig4, ax4 = plt.subplots(1, 1)
ax4.scatter(atoms, proj_e1, linewidth=1.0) 
ax4.set_xlabel("Atom index", size="medium")
ax4.set_ylabel("Atomic projection", size="medium")
ax4.set_title("VI(0) 8d-I4c, saddle, e1 (acoustic mode)")
fig4.tight_layout()
fig4.savefig("./atom_proj_e1.png")

fig5, ax5 = plt.subplots(1, 1)
ax5.scatter(atoms, proj_e4, linewidth=1.0) 
ax5.set_xlabel("Atom index", size="medium")
ax5.set_ylabel("Atomic projection", size="medium")
ax5.set_title("VI(0) 8d-I4c, saddle, e1 (acoustic mode)")
fig5.tight_layout()
fig5.savefig("./atom_proj_e4.png")