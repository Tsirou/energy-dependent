import os, sys

import numpy as np

files      = [ "SED_ECPL_pts.txt", "SED_PL_pts.txt", "SED_STD_ECPL_pts.txt", "SED_STD_PL_pts.txt"]
path       = "/home/tsirou/Documents/Analyses/MSH_15-52/energy-dependent/MultiEnergyBins/spectra/"

for filename in files:

    energy, flux, flux_min, flux_max = np.loadtxt(path + filename, dtype=float, unpack=True)

    ecsv_file   = open(path + filename[:-3] + "ecsv", "w")
    ecsv_file.write("# %ECSV 0.9\n\
# ---\n\
# datatype:\n\
# - {name: energy, unit: TeV, datatype: float64}\n\
# - {name: flux, unit: 1 / (TeV cm2 s), datatype: float64}\n\
# - {name: flux_error_lo, unit: 1 / (TeV cm2 s), datatype: float64}\n\
# - {name: flux_error_hi, unit: 1 / (TeV cm2 s), datatype: float64}\n\
energy flux flux_error_lo flux_error_hi\n")

    for pt in range(0, len(energy)):
        ecsv_file.write(str(energy[pt]) + " " + \
                        str(flux[pt]) + " " + \
                        str(flux_min[pt]) + " " + \
                        str(flux_max[pt]) + "\n")

#
# files      = [ "spectrum_2005.txt"]
# path       = "/home/tsirou/Documents/Analyses/MSH_15-52/energy-dependent/MultiEnergyBins/spectra/"
#
# for filename in files:
#
#     energy_min, energy_max, energy_mean, flux, flux_error = np.loadtxt(path + filename, dtype=float, unpack=True, skiprows=1)
#
#     ecsv_file   = open(path + filename[:-3] + "ecsv", "w")
#     ecsv_file.write("# %ECSV 0.9\n\
# # ---\n\
# # datatype:\n\
# # - {name: energy, unit: TeV, datatype: float64}\n\
# # - {name: flux, unit: 1 / (TeV cm2 s), datatype: float64}\n\
# # - {name: flux_error_lo, unit: 1 / (TeV cm2 s), datatype: float64}\n\
# # - {name: flux_error_hi, unit: 1 / (TeV cm2 s), datatype: float64}\n\
# energy flux flux_error_lo flux_error_hi\n")
#
#     for pt in range(0, len(energy_mean)):
#         ecsv_file.write(str(energy_mean[pt]) + " " + \
#                         str(flux[pt]) + " " + \
#                         str(flux[pt] - flux_error[pt]) + " " + \
#                         str(flux[pt] + flux_error[pt]) + "\n")

