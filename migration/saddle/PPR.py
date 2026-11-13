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
        #sum = sum + mass[k]*mod_sq(disp[k])
        sum = sum + mod_sq(disp[k])
    return sum*sum


# for the 2nd energy moment, use squares of the inputs
def M2(mass, disp):
    N_atoms = len(disp)
    sum = 0
    for k in range(N_atoms):
        #sum = sum + (mass[k]*mod_sq(disp[k]))**2
        sum = sum + (mod_sq(disp[k]))*(mod_sq(disp[k]))
    return sum  

## calculate PPR ##
m = get_mass(data_atoms)

# very first mode is the unstable negative freq, highly localised on migrating ion
# wants to move from local max of PES (TS) to local min 
e0 = get_evec(data_bands, 0)

# acoustic modes [0,2] aren't relevant; at q=0 they should all be collective translations
# so PPR ~ 1
e1 = get_evec(data_bands, 1)

# first optical mode, probably more localised? so PPR < 1
e4 = get_evec(data_bands, 4)

p0 = M1(m, e0)/(len(e0)*M2(m, e0))
p1 = M1(m, e1)/(len(e1)*M2(m, e1))
p4 = M1(m, e4)/(len(e4)*M2(m, e4))

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
        #print("i = ", i)
        e0_full[i] = e0[k][j]
        #print(e0_full[i])
        count = count + 1
#print(e0_full)

print("|e0|^2 = ", mod_sq(e0_full))

# loop over all phonon modes j, plot the PPRs?
# finally make an Natoms variable, won't be too slow here because not referencing in different scope
Nj = len(data_bands)
Nk = len(m)
freqs = [None]*Nj
modes = [None]*Nj
pj = [None]*Nj
for j in range(Nj):
    ej = get_evec(data_bands, j)
    pj[j] = M1(m, ej)/(Nk*M2(m, ej)) 
    freqs[j] = data_bands[j]["frequency"] 
    modes[j] = j+1
    
# maybe plot a delta function at each frequency instead??
# or plot j instead of freq for even spacing
    
fig1, ax1 = plt.subplots(1, 1)
#ax1.scatter(freqs, pj, marker='o', linewidths=0.1) 
ax1.scatter(freqs, pj, linewidth=1.0) 
#ax1.set_xlim()
#ax1.set_ylim()
ax1.set_xlabel("Phonon frequency (THz)", size="medium")
ax1.set_ylabel("Phonon participation ratio", size="medium")
#ax1.legend()
ax1.set_title("VI(0) 8d-I4c, saddle")
fig1.tight_layout()
fig1.savefig("./PPR_freq.png")

fig2, ax2 = plt.subplots(1, 1)
#ax1.scatter(freqs, pj, marker='o', linewidths=0.1) 
ax2.scatter(modes, pj, linewidth=1.0) 
#ax1.set_xlim()
#ax1.set_ylim()
ax2.set_xlabel("Phonon mode index j", size="medium")
ax2.set_ylabel("Phonon participation ratio", size="medium")
#ax1.legend()
ax2.set_title("VI(0) 8d-I4c, saddle")
fig2.tight_layout()
fig2.savefig("./PPR_j.png")