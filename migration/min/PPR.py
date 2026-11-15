import numpy as np
import yaml 
from yaml import CLoader as Loader, CDumper as Dumper
import matplotlib.pyplot as plt

# READ FROM BANDS...NOT QPOINTS...NEED MASS DATA!
with open("./band.yaml", 'r') as f:
#    data = yaml.load(f, Loader=yaml.SafeLoader)
    data = yaml.load(f, Loader=Loader)

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

#m = get_mass(data_atoms)

## calculate PPR ##

# acoustic modes [0,2] aren't relevant; at q=0 they should all be collective translations
# so PPR ~ 1
e0 = get_evec(data_bands, 0)

# first optical mode, probably more localised? so PPR < 1
e3 = get_evec(data_bands, 3)

p0 = M1(e0)/(len(e0)*M2(e0))
ip0 = 1/M2(e0)
p3 = M1(e3)/(len(e3)*M2(e3))
ip3 = 1/M2(e3)

print("(sum of |n|^2)^2 (acoustic) = ", M1(e0))
print("sum of (|n|^2)^2 (acoustic) = ", M2(e0))

print("PPR for band 1 (acoustic)", "at q = 0:", p0)
#print("IPPR for band 1 (acoustic)", "at q = 0:", ip0)
print("PPR for band 4 (optical)", "at q = 0:", p3)
#print("IPPR for band 2 (acoustic)", "at q = 0:", ip1)

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
proj_Pb = [None]*Nj
proj_Cs = [None]*Nj
proj_I = [None]*Nj
proj_Br = [None]*Nj
#ipj = [None]*Nj
for j in range(Nj):
    ej = get_evec(data_bands, j)
    pj[j] = M1(ej)/(Nk*M2(ej)) 
    #ipj[j] = 1/M2(ej) 
    freqs[j] = data_bands[j]["frequency"] 
    modes[j] = j+1
    
    
# also plot the magnitude of atom projections, for each mode
# each atomic displacement component r_k
# plot the projection into each of 159 atoms for all 477 modes??
# or show proj into atomic species Cs, Pb, I, Br
# since I only care about behaviour of individual modes, can look at all atoms
proj_e0 = [None]*Nk
proj_e3 = [None]*Nk
atoms = [None]*Nk
for k in range(Nk):
    proj_e0[k] = mod_sq(e0[k])
    proj_e3[k] = mod_sq(e3[k])
    atoms[k] = k+1
    
# maybe plot a delta function at each frequency instead??
# or plot j instead of freq for even spacing
    
fig1, ax1 = plt.subplots(1, 1)
ax1.scatter(freqs, pj, linewidth=1.0) 
ax1.set_xlabel("Phonon frequency (THz)", size="medium")
ax1.set_ylabel("Phonon participation ratio", size="medium")
ax1.set_title("VI(0) 8d-I4c, min")
fig1.tight_layout()
fig1.savefig("./PPR_freq.png")

fig2, ax2 = plt.subplots(1, 1)
ax2.scatter(modes, pj, linewidth=1.0) 
ax2.set_xlabel("Phonon mode index j", size="medium")
ax2.set_ylabel("Phonon participation ratio", size="medium")
ax2.set_title("VI(0) 8d-I4c, min")
fig2.tight_layout()
fig2.savefig("./PPR_j.png")

fig3, ax3 = plt.subplots(1, 1)
ax3.scatter(atoms, proj_e0, linewidth=1.0) 
ax3.set_xlabel("Atom index", size="medium")
ax3.set_ylabel("Atomic projection", size="medium")
ax3.set_title("VI(0) 8d-I4c, min, e0 (acoustic mode)")
fig3.tight_layout()
fig3.savefig("./atom_proj_e0.png")

fig4, ax4 = plt.subplots(1, 1)
ax4.scatter(atoms, proj_e3, linewidth=1.0) 
ax4.set_xlabel("Atom index", size="medium")
ax4.set_ylabel("Atomic projection", size="medium")
ax4.set_title("VI(0) 8d-I4c, min, e3 (1st optical mode)")
fig4.tight_layout()
fig4.savefig("./atom_proj_e3.png")