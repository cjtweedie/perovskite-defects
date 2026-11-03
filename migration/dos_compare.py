import numpy as np
import matplotlib.pyplot as plt

data_min = np.loadtxt("./total_dos_min.dat")
data_sad = np.loadtxt("./total_dos_saddle.dat")

freq_min = data_min[:,0]
dos_min = data_min[:,1]
freq_sad = data_sad[:,0]
dos_sad = data_sad[:,1]


fig1, ax1 = plt.subplots(1, 1)
ax1.plot(freq_min, dos_min, label="Min")
ax1.plot(freq_sad, dos_sad, label="Saddle")
ax1.axhline(0, c="black", linestyle='--', linewidth=0.9)    
ax1.set_xlim(freq_sad[0], max(freq_min[-1], freq_sad[-1]))
ax1.set_ylim(0, max(np.max(dos_min), np.max(dos_sad))+5)
#ax1.set_xticks(kmag_ticks)
#ax1.set_xticklabels(kpath, size="medium")
#ax1.tick_params(axis="both", direction="in")
ax1.set_xlabel("Frequency [THz]", size="medium")
ax1.set_ylabel("DOS [a.u.]", size="medium")
ax1.legend()
fig1.tight_layout()
fig1.savefig("./dos_compare.png")