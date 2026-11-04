import numpy as np
import matplotlib.pyplot as plt

data_V0_min = np.loadtxt("./VI0/Pb57/8d-I4c/min/total_dos.dat")
data_V0_sad = np.loadtxt("./VI0/Pb57/8d-I4c/saddle/total_dos.dat")
data_Vp_min = np.loadtxt("./VI+/Pb57/8d-I4c/min/total_dos.dat")
data_Vp_sad = np.loadtxt("./VI+/Pb57/8d-I4c/saddle/total_dos.dat")

freq_V0_min = data_V0_min[:,0]
dos_V0_min = data_V0_min[:,1]
freq_V0_sad = data_V0_sad[:,0]
dos_V0_sad = data_V0_sad[:,1]

freq_Vp_min = data_Vp_min[:,0]
dos_Vp_min = data_Vp_min[:,1]
freq_Vp_sad = data_Vp_sad[:,0]
dos_Vp_sad = data_Vp_sad[:,1]

dos_V0_diff = np.zeros(len(freq_V0_min))
dos_Vp_diff = np.zeros(len(freq_Vp_min))
dos_V_diff = np.zeros(len(freq_V0_min))
f0 = open("total_dos_V0_diff.txt", 'w')
fp = open("total_dos_Vp_diff.txt", 'w')
start_diff = np.where(freq_V0_sad == freq_V0_min[0])[0][0]
for j in range(len(freq_V0_min)):
    dos_V0_diff[j] = dos_V0_min[j]-dos_V0_sad[j+start_diff]
    dos_Vp_diff[j] = dos_Vp_min[j]-dos_Vp_sad[j+start_diff]
    dos_V_diff[j] = dos_V0_diff[j] - dos_Vp_diff[j]
    print(freq_V0_min[j], dos_V0_diff[j], file=f0)
    print(freq_Vp_min[j], dos_Vp_diff[j], file=fp)


fig1, ax1 = plt.subplots(1, 1)
ax1.plot(freq_V0_min, dos_V0_diff, label="VI0")   
ax1.plot(freq_Vp_min, dos_Vp_diff, label="VI+") 
ax1.plot(freq_V0_min, 0.5*dos_V_diff, label="VI0 - VI+")
ax1.plot(freq_V0_min, np.zeros(len(freq_V0_min)), color="black", linestyle='--')
ax1.set_xlim(freq_V0_min[0], max(freq_V0_min[-1], freq_V0_sad[-1]))
ax1.set_ylim(min(np.min(dos_V0_diff), np.min(dos_Vp_diff)), max(np.max(dos_V0_diff), np.max(dos_Vp_diff)))
ax1.set_xlabel("Frequency [THz]", size="medium")
ax1.set_ylabel("Density of states [a.u.]", size="medium")
ax1.legend()
ax1.set_title("Phonon DOS diff for VI+, VI0")
fig1.tight_layout()
fig1.savefig("./dos_compare_charge.png")