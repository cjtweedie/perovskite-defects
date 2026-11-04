import numpy as np
import matplotlib.pyplot as plt

data_min = np.loadtxt("./min/total_dos.dat")
data_sad = np.loadtxt("./saddle/total_dos.dat")

freq_min = data_min[:,0]
dos_min = data_min[:,1]
freq_sad = data_sad[:,0]
dos_sad = data_sad[:,1]

dos_diff = np.zeros(len(freq_min))
f = open("total_dos_diff.txt", 'w')
start_diff = np.where(freq_sad == freq_min[0])[0][0]
for j in range(len(freq_min)):
    dos_diff[j] = dos_min[j]-dos_sad[j+start_diff]
    print(freq_min[j], dos_diff[j], file=f)


fig1, ax1 = plt.subplots(1, 1)
ax1.plot(freq_min, dos_min, label="min")
ax1.plot(freq_sad, dos_sad, label="saddle")
ax1.plot(freq_min, 10*dos_diff, label="a*(min - saddle)")   
ax1.plot(freq_min, np.zeros(len(freq_min)), color="black", linestyle='--')
ax1.set_xlim(freq_min[0], max(freq_min[-1], freq_sad[-1]))
ax1.set_ylim(10*np.min(dos_diff)-3, max(np.max(dos_min), np.max(dos_sad))+3)
ax1.set_xlabel("Frequency [THz]", size="medium")
ax1.set_ylabel("Density of states [a.u.]", size="medium")
ax1.legend()
ax1.set_title("Min vs saddle phonon DOS at q=0")
fig1.tight_layout()
fig1.savefig("./dos_compare.png")