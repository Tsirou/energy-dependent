import os
import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from heapq import merge
from display import slice_energy_plot, slice_label_excess_plot, reading_energy_bins

B_bins, B_on, B_off                = np.loadtxt(names.path_multiEbins + names.filename_Bsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
B_E_min, B_E_max, B_alpha, B_gamma = np.loadtxt(names.path_multiEbins + names.filename_Bsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

S_bins, S_on, S_off                = np.loadtxt(names.path_multiEbins + names.filename_Ssteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
S_E_min, S_E_max, S_alpha, S_gamma = np.loadtxt(names.path_multiEbins + names.filename_Ssteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

P_bins, P_on, P_off                = np.loadtxt(names.path_multiEbins + names.filename_Psteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
P_E_min, P_E_max, P_alpha, P_gamma = np.loadtxt(names.path_multiEbins + names.filename_Psteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

N_bins, N_on, N_off                = np.loadtxt(names.path_multiEbins + names.filename_Nsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
N_E_min, N_E_max, N_alpha, N_gamma = np.loadtxt(names.path_multiEbins + names.filename_Nsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

H_bins, H_on, H_off                = np.loadtxt(names.path_multiEbins + names.filename_Hsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
H_E_min, H_E_max, H_alpha, H_gamma = np.loadtxt(names.path_multiEbins + names.filename_Hsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

NH_bins, NH_on, NH_off                = np.loadtxt(names.path_multiEbins + names.filename_NHsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
NH_E_min, NH_E_max, NH_alpha, NH_gamma = np.loadtxt(names.path_multiEbins + names.filename_NHsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

Hn_bins, Hn_on, Hn_off                = np.loadtxt(names.path_multiEbins + names.filename_Hnsteps, dtype=int, skiprows=1, usecols=(0,3,4), unpack=True)
Hn_E_min, Hn_E_max, Hn_alpha, Hn_gamma = np.loadtxt(names.path_multiEbins + names.filename_Hnsteps, dtype=float, skiprows=1, usecols=(1,2,5,6), unpack=True)

colors                             = ["cornflowerblue", "springgreen", "orange", "indianred","crimson","indigo","navy","teal","limegreen","brown","magenta","darkgrey"]



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


# # plt.plot( S_E_min, S_on, color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ [S]')
# # plt.plot( S_E_max, S_on, color='darkgreen', marker=".", drawstyle='steps-pre', linestyle='solid')
# #
# #
# # plt.plot( S_E_min, S_off / S_alpha , color='green', marker=".", linestyle='solid', drawstyle='steps-post', label='a$^{-1}$ N$_{OFF}$ [S]')
# # plt.plot( S_E_max, S_off / S_alpha , color='green', marker=".", drawstyle='steps-pre', linestyle='solid')
# #
# # plt.plot( S_E_min, S_off , color='limegreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{OFF}$ [S]')
# # plt.plot( S_E_max, S_off , color='limegreen', marker=".", drawstyle='steps-pre', linestyle='solid')
# #
#
# plt.plot( S_E_min, S_on / (S_off * S_alpha) , color='green', marker=".", linestyle='solid', drawstyle='steps-post', label='a$^{-1}$ N$_{OFF}$ [S]')
# plt.plot( S_E_max, S_on / (S_off * S_alpha) , color='green', marker=".", drawstyle='steps-pre', linestyle='solid')
#
# plt.legend(fontsize=10,numpoints=1)
#
# plt.show()
#
# os._exit(0)



print("\nHiRes previous steps bin choice")

nb_bins           = 4
H_computed_slices_E, H_computed_slices_ex = slicing_energy_bins(H_on, H_off, H_alpha, H_E_min, H_E_max, nb_bins)

plt.title('Energy vs excess\n\n Bins : [0.3 - 0.6 || 0.6 - 0.9 || 0.9 - 3.0 || 3.0 - 100] TeV')
plt.xlabel('Energy (TeV)')
plt.ylabel('Excess')

P_err_exc   = np.sqrt(P_on) + (np.sqrt(P_off) / P_alpha )
H_err_exc  = np.sqrt(H_on) - (np.sqrt(H_off) / H_alpha )

plt.semilogx( P_E_min, P_on - (P_off / P_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [std]')
plt.semilogx( P_E_max, P_on - (P_off / P_alpha ), color='darkgreen', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( P_E_min, (P_on - (P_off / P_alpha ) ) + P_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( P_E_max, (P_on - (P_off / P_alpha ) ) + P_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( P_E_min, (P_on - (P_off / P_alpha ) ) - P_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( P_E_max, (P_on - (P_off / P_alpha ) ) - P_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)


plt.semilogx( H_E_min, H_on - (H_off / H_alpha ), color='darkorange', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [HiRes]')
plt.semilogx( H_E_max, H_on - (H_off / H_alpha ), color='darkorange', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( H_E_min, (H_on - (H_off / H_alpha ) ) + H_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( H_E_max, (H_on - (H_off / H_alpha ) ) + H_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( H_E_min, (H_on - (H_off / H_alpha ) ) - H_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( H_E_max, (H_on - (H_off / H_alpha ) ) - H_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)


E_cuts  = [0.3, 0.6, 0.9, 3.0, 100.]

exc     = 0.
Ex_cuts = []
Ex_cuts.append(0.)
for i in range(0, len(N_on)):
    exc = exc + (P_on[i] - P_off[i] / P_alpha[i])
    Ex_cuts.append(exc)

slice_energy_plot(E_cuts, colors)
slice_label_excess_plot(E_cuts, Ex_cuts, P_gamma, colors)


exc     = 0.
Ex_cuts = []
Ex_cuts.append(0.)
for i in range(0, len(H_on)):
    exc = exc + (H_on[i] - H_off[i] / H_alpha[i])
    Ex_cuts.append(exc)

slice_energy_plot(E_cuts, colors)
slice_label_excess_plot(E_cuts, Ex_cuts, H_gamma, colors)

plt.legend(fontsize=10,numpoints=1,loc=3)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_s'  + '.png', dpi=1000)


plt.clf()

print("\nNew steps bin choice")

nb_bins           = 4
N_computed_slices_E, N_computed_slices_ex = slicing_energy_bins(N_on, N_off, N_alpha, N_E_min, N_E_max, nb_bins)

plt.title('Energy vs excess\n\n Bins : [0.3 - 0.6 || 0.6 - 1.2 || 1.2 - 2.4 || 2.4 - 100] TeV')
plt.xlabel('Energy (TeV)')
plt.ylabel('Excess')


N_err_exc   = np.sqrt(N_on) + (np.sqrt(N_off) / N_alpha )
Hn_err_exc  = np.sqrt(Hn_on) - (np.sqrt(Hn_off) / Hn_alpha )

plt.semilogx( N_E_min, N_on - (N_off / N_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [std]')
plt.semilogx( N_E_max, N_on - (N_off / N_alpha ), color='darkgreen', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( N_E_min, (N_on - (N_off / N_alpha ) ) + N_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( N_E_max, (N_on - (N_off / N_alpha ) ) + N_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( N_E_min, (N_on - (N_off / N_alpha ) ) - N_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( N_E_max, (N_on - (N_off / N_alpha ) ) - N_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)


plt.semilogx( Hn_E_min, Hn_on - (Hn_off / Hn_alpha ), color='darkorange', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [HiRes]')
plt.semilogx( Hn_E_max, Hn_on - (Hn_off / Hn_alpha ), color='darkorange', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( Hn_E_min, (Hn_on - (Hn_off / Hn_alpha ) ) + Hn_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( Hn_E_max, (Hn_on - (Hn_off / Hn_alpha ) ) + Hn_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( Hn_E_min, (Hn_on - (Hn_off / Hn_alpha ) ) - Hn_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( Hn_E_max, (Hn_on - (Hn_off / Hn_alpha ) ) - Hn_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)


E_cuts  = [0.3, 0.6, 1.2, 2.4, 100.]

exc     = 0.
Ex_cuts = []
Ex_cuts.append(0.)
for i in range(0, len(N_on)):
    exc = exc + (N_on[i] - N_off[i] / N_alpha[i])
    Ex_cuts.append(exc)

slice_energy_plot(E_cuts, colors)
slice_label_excess_plot(E_cuts, Ex_cuts, N_gamma, colors)


exc     = 0.
Ex_cuts = []
Ex_cuts.append(0.)
for i in range(0, len(Hn_on)):
    exc = exc + (Hn_on[i] - Hn_off[i] / Hn_alpha[i])
    Ex_cuts.append(exc)

slice_energy_plot(E_cuts, colors)
slice_label_excess_plot(E_cuts, Ex_cuts, Hn_gamma, colors)

plt.legend(fontsize=10,numpoints=1,loc=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_s'  + '.png', dpi=1000)

plt.clf()



print("\nSmall log steps bin choice")

plt.title('Energy vs excess')
plt.xlabel('Energy (TeV)')
plt.ylabel('Excess')


S_err_exc   = np.sqrt(S_on) + (np.sqrt(S_off) / S_alpha )
NH_err_exc  = np.sqrt(NH_on) - (np.sqrt(NH_off) / NH_alpha )

plt.semilogx( S_E_min, S_on - (S_off / S_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [std]')
plt.semilogx( S_E_max, S_on - (S_off / S_alpha ), color='darkgreen', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( S_E_min, (S_on - (S_off / S_alpha ) ) + S_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( S_E_max, (S_on - (S_off / S_alpha ) ) + S_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( S_E_min, (S_on - (S_off / S_alpha ) ) - S_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( S_E_max, (S_on - (S_off / S_alpha ) ) - S_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)


plt.semilogx( NH_E_min, NH_on - (NH_off / NH_alpha ), color='darkorange', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a$^{-1}$ N$_{OFF}$ [HiRes]')
plt.semilogx( NH_E_max, NH_on - (NH_off / NH_alpha ), color='darkorange', marker=".", drawstyle='steps-pre', linestyle='solid')

plt.semilogx( NH_E_min, (NH_on - (NH_off / NH_alpha )) + NH_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( NH_E_max, (NH_on - (NH_off / NH_alpha )) + NH_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( NH_E_min, (NH_on - (NH_off / NH_alpha )) - NH_err_exc, color='darkorange', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( NH_E_max, (NH_on - (NH_off / NH_alpha )) - NH_err_exc, color='darkorange', drawstyle='steps-pre', linestyle='--', alpha = 0.5)

reading_energy_bins(S_on, S_off, S_alpha, S_E_min, S_E_max, colors)
reading_energy_bins(NH_on, NH_off, NH_alpha, NH_E_min, NH_E_max, colors)

plt.legend(fontsize=10,numpoints=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_s'  + '.png', dpi=1000)

plt.clf()


print("\nBig log steps bin choice")

plt.title('Energy vs excess\n\n Log steps [0.3 - 75.4] TeV')
plt.xlabel('Energy (TeV)')
plt.ylabel('Excess')

B_err_exc   = np.sqrt(B_on) + (np.sqrt(B_off) / B_alpha )

plt.semilogx( B_E_min , B_on - (B_off / B_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-post', label='N$_{ON}$ - a $^{-1}$ N$_{OFF}$ [std]')
plt.semilogx( B_E_max , B_on - (B_off / B_alpha ), color='darkgreen', marker=".", linestyle='solid', drawstyle='steps-pre')

plt.semilogx( B_E_min, (B_on - (B_off / B_alpha ) ) + B_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( B_E_max, (B_on - (B_off / B_alpha ) ) + B_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)
plt.semilogx( B_E_min, (B_on - (B_off / B_alpha ) ) - B_err_exc, color='darkgreen', linestyle='--', drawstyle='steps-post', alpha = 0.5)
plt.semilogx( B_E_max, (B_on - (B_off / B_alpha ) ) - B_err_exc, color='darkgreen', drawstyle='steps-pre', linestyle='--', alpha = 0.5)

reading_energy_bins(B_on, B_off, B_alpha, B_E_min, B_E_max, colors)


plt.legend(fontsize=10,numpoints=1)
plt.show()
#plt.savefig(names.path_multiEbins_results + 'energy_VS_excess_b'  + '.png', dpi=1000)

plt.clf()




