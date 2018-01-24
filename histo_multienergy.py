import os
import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from heapq import merge
from display import slice_energy_plot, slice_label_excess_plot

B_bins, B_on, B_off                = np.loadtxt(names.path_multiEbins + names.filename_Bsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
B_E_min, B_E_max, B_alpha, B_gamma = np.loadtxt(names.path_multiEbins + names.filename_Bsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

S_bins, S_on, S_off                = np.loadtxt(names.path_multiEbins + names.filename_Ssteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
S_E_min, S_E_max, S_alpha, S_gamma = np.loadtxt(names.path_multiEbins + names.filename_Ssteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

P_bins, P_on, P_off                = np.loadtxt(names.path_multiEbins + names.filename_Psteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
P_E_min, P_E_max, P_alpha, P_gamma = np.loadtxt(names.path_multiEbins + names.filename_Psteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

colors                             = ["cornflowerblue", "springgreen", "orange", "indianred","crimson","indigo","navy","teal","limegreen","yellow","brown","silver"]



def slicing_energy_bins(X_on, X_off, X_alpha, X_E_min, X_E_max, X_nb_bins):

    excess_sum = sum((X_on - (X_off / X_alpha)))

    slices = []

    for j in range(0, X_nb_bins):
        slices.append(excess_sum / X_nb_bins * j)

    exc                  = 0.
    excess               = []
    X_computed_slices_E  = []
    X_computed_slices_ex = []

    for i in range(0, len(X_on)):
        exc = exc + (X_on[i] - (X_off[i] / X_alpha[i]))
        excess.append(exc)

    X_computed_slices_E.append(X_E_min[0])
    X_computed_slices_ex.append(0.)
    for i in range(0, len(X_on) - 1):
        for s in range(0, len(slices) - 1):
            if (excess[i] >= slices[s] and excess[i] < slices[s + 1] and excess[i + 1] > slices[s + 1]):
                print "SLICE IT UP : ", excess[i], slices[s + 1], X_E_max[i]
                X_computed_slices_E.append(X_E_max[i])
                X_computed_slices_ex.append(excess[i])

    X_computed_slices_E.append(X_E_max[-1])
    X_computed_slices_ex.append(excess_sum)
    return X_computed_slices_E, X_computed_slices_ex


###############################################################################################################
###############################################################################################################

plt.clf()

print("\nSmall log steps bin choice")

nb_bins           = 12
S_computed_slices_E, S_computed_slices_ex = slicing_energy_bins(S_on, S_off, S_alpha, S_E_min, S_E_max, nb_bins)

plt.title('Energy vs excess')
plt.xlabel('Energy (TeV)')
plt.ylabel('Excess')

plt.xticks(np.arange(B_E_min[0], B_E_max[-1], 0.3))

log_bins = list(merge(B_E_min,B_E_max))[::2]

plt.semilogx( S_E_min, S_on - (S_off / S_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [S]')
plt.semilogx( S_E_max, S_on - (S_off / S_alpha ), color='darkgreen', marker=".", drawstyle='steps-pre', linestyle='solid')

slice_energy_plot(E_cuts=S_computed_slices_E, colors=colors)
slice_label_excess_plot(S_computed_slices_E, S_computed_slices_ex, S_gamma, colors)


plt.legend(fontsize=10,numpoints=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_s'  + '.png', dpi=1000)

plt.clf()

print("\nBig log steps bin choice")

nb_bins           = 3
B_computed_slices_E, B_computed_slices_ex = slicing_energy_bins(B_on, B_off, B_alpha, B_E_min, B_E_max, nb_bins)

plt.semilogx( B_E_min , B_on - (B_off / B_alpha ), color='darkblue', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a $^{-1}$ N$_{OFF}$ [B]')
plt.semilogx( B_E_max , B_on - (B_off / B_alpha ), color='darkblue', marker=".", linestyle='solid', drawstyle='steps-pre')

slice_energy_plot(E_cuts=B_computed_slices_E, colors=colors)
slice_label_excess_plot(B_computed_slices_E, B_computed_slices_ex, B_gamma, colors)


plt.legend(fontsize=10,numpoints=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_b'  + '.png', dpi=1000)

plt.clf()


print("\nPrevious bin choice")

nb_bins    = 4
P_computed_slices_E, P_computed_slices_ex = slicing_energy_bins(P_on, P_off, P_alpha, P_E_min, P_E_max, nb_bins)


plt.semilogx( P_E_min , P_on - (P_off / P_alpha ), color='darkred', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a $^{-1}$ N$_{OFF}$ [P]')
plt.semilogx( P_E_max , P_on - (P_off / P_alpha ), color='darkred', marker=".", linestyle='solid', drawstyle='steps-pre')


E_cuts  = [0.3, 0.6, 0.9, 3.0, 100.]
exc     = 0.
Ex_cuts = []
Ex_cuts.append(0.)
for i in range(0, len(P_on)):
    exc = exc + (P_on[i] - (P_off[i] / P_alpha[i]))
    Ex_cuts.append(exc)
#Ex_cuts.append(sum((P_on - (P_off / P_alpha))))

slice_energy_plot(E_cuts, colors)
slice_label_excess_plot(E_cuts, Ex_cuts, P_gamma, colors)


plt.legend(fontsize=10,numpoints=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_p'  + '.png', dpi=1000)


