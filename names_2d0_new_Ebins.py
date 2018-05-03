import sys, os
import numpy as np


fit_psf           = ""
rebinned          = ""
psf               = ""
confidence        = ""
psf_convl         = ""
save_option       = ""
fit_option        = ""
bkg_option        = ""

#==============================================================================
# Input arguments :
# 0               : python [command]
# 1               : analysis [names]
# 2               : component_list_index_start [main]
# 3               : component_list_index_end [main]
# 4               : sliced_selection [names]
# 5               : nb_of_slice [names]
#==============================================================================

#==============================================================================
# Various variables
#==============================================================================


random_convert_shift =  0  # 1 or 0 --- Depends on the matrix conventions
display_ds9          = 'n' # yes (y) or no (n) for displaying with ds9 the fit steps
circular_selection   = 'n' # yes (y) or no (n) for selecting a circular region for fitting
sliced_selection     = sys.argv[4] # yes (y) or no (n) for selecting a circular region for fitting
bkg_thaw             = 'n' # yes (y) or no (n) for thawing the background amplitude
bkg_amplitude        = 1.0 # usually set to 1.0 but could also be +/- 2%
confidence_levels    = 'n' # yes (y) or no (n) for computing the confidence intervals (time-consuming!)


r_psf_0           = 0
r_psf_1           = 0
r_psf_2           = 0


r_obs             = 0
r_intrinsic       = 0

radius_it         = 0

all_sum           = 0.
sum_68            = 0.

mu                = 0.
sigma             = 0.

factor_pix2deg    = 0.


###############################################
analysis          = int(sys.argv[1])
# Tag for the different analyses :
#     0    : PA 2.0deg 0.3 - 0.6 TeV
#     1    : PA 2.0deg 0.6 - 1.2 TeV
#     2    : PA 2.0deg 1.2 - 2.4 TeV
#     3    : PA 2.0deg > 2.4 TeV
#     4    : PA 2.0deg Std
###############################################

nb_pa_analyses    = 5

pulsar_then_rebin = 'n' # Check in case one should multiply the x-ray templates before rebinning them --Not fully implemented || Also used for the
gif_norm          = 'n' # For the normalisation of the colorbar in the gifs generation

psf_fudge         = -1  # 0.02 or -1 value in sigma degrees! Has to be converted into pixels. If there is no fudge factor then insert -1

method            = "simplex" # "moncar"

dim_for_a         = '400x400'
# For the hdu FITS extensions for the gamma-like maps
hdu_ext = 1  # 0 for gammapy maps, 1 for PA

# psf_fudge and bkg_amplitude arguments for the systematics estimation


# Display GIF normalisation for the colorbar
cb_extrema              = np.ndarray((nb_pa_analyses, 2))
cb_extrema[0,:]         = [1.0, 1.8]
cb_extrema[1,:]         = [1.0, 1.6]
cb_extrema[2,:]         = [1.0, 1.6]
cb_extrema[3,:]         = [1.0, 1.2]
cb_extrema[4,:]         = [1.0, 4.0]


cb_resids               = np.ndarray((nb_pa_analyses, 2))
cb_resids[0,:]          = [-3.0, 3.0]
cb_resids[1,:]          = [-3.0, 3.0]
cb_resids[2,:]          = [-3.0, 3.0]
cb_resids[3,:]          = [-3.0, 3.0]
cb_resids[4,:]          = [-4.0, 4.0]


#For the exposure_gamma.py script
analysis_energy_low        = [0.3, 0.6, 1.2, 2.4, 0.3]
analysis_energy_high       = [0.6, 1.2, 2.4, 100., 100.]


source            = ["MSH15-52---0.3_0.6TeV", "MSH15-52---0.6_1.2TeV","MSH15-52---1.2-2.4TeV", "MSH15-52---gt2.4TeV", "MSH15-52---std"]
tags              = ["0-3_0-6TeV", "0-6_1-2TeV", "1-2_2-4TeV", "gt2-4TeV", "std"]
#fit_component     = ["xr0-ScP","xr0-Sc", "xr-ScP", "xr-Sc"]
fit_component     = ["xr0","xr0-2G","xr","xr-2G","xr-G1","xr-sh","xr-dc","xr0-ScP","xr0-Sc", "xr-ScP", "xr-Sc"]
#fit_component     = ["xr_psf","xr_psf","xr-2G_psf","xr-G1_psf","xr-sh_psf","xr-dc_psf"]
#fit_component     = ["xr-2G_bkg_d"]



alpha_factor_pa   = [0.0, 0.0, 1.3, 1.3, 1.3, 1.3, 1.3, 0.0, 0.0, 1.3, 1.3] # when alpha is frozen to ~ the best fit value for the Std cut model + 2G
#alpha_factor_pa   = [0.0, 0.0, 1.3, 1.3] # for the Sersic profile investigation



directory         = ["energy-dependent/MultiEnergyBins/new/Ebins/"] * nb_pa_analyses

path_data         = "/home/tsirou/Documents/Analyses/DATA/"
path_g            = ["/home/tsirou/Documents/Analyses/MSH_15-52/"] * (nb_pa_analyses)
path              = [""] * (nb_pa_analyses)
for n in range(0,nb_pa_analyses):
    path[n]   = path_g[n] + directory[n]

path_x            = "/home/tsirou/Documents/Analyses/Xray_MSH_15-52/"

path_psf          = [p + "psf/" for p in path]
path_template     = [p + "templates/" for p in path]
path_expo         = [p + "exposure/" for p in path]
path_runs         = [p + "runs/" for p in path]

save_path         = [p + "Fit/" for p in path]
save_path_gif     = [s + 'GIF/' for s in save_path]

if (sliced_selection.find('y') != -1):
    save_path_sliced  = [p + "sliced/" for p in path]
    save_path         = [p + "sliced/Fit/" for p in path]
    save_path_gif     = [s + 'GIF/' for s in save_path]


# Path for all the results txt files (global)
results_path      = "/home/tsirou/Documents/Analyses/Results/Energy-dependent/MultiEnergyBins/2d0_new/"
if (sliced_selection.find('y') != -1):
    results_path   = results_path + "sliced/"


filename_runs     = ["Results_MSH_15_52_0_3_0_6_TeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_0_6_1_2_TeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_1_2_2_4_TeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_gt2_4_TeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_std_2d0_ModelPlus_Std_mergedrun_list.txt"]

filename_bkg      = ["Results_MSH_15_52_0_3_0_6_TeV_Custom_mergedbackground.fits", "Results_MSH_15_52_0_6_1_2_TeV_Custom_mergedbackground.fits", "Results_MSH_15_52_1_2_2_4_TeV_Custom_mergedbackground.fits", "Results_MSH_15_52_1_2_2_4_TeV_Custom_mergedbackground.fits", "Results_MSH_15_52_std_2d0_ModelPlus_Std_mergedbackground.fits"]
filename_gamma    = ["Results_MSH_15_52_0_3_0_6_TeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_0_6_1_2_TeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_1_2_2_4_TeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_gt2_4_TeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_std_2d0_ModelPlus_Std_mergedgamma-like.fits"]
filename_excess   = ["Results_MSH_15_52_0_3_0_6_TeV_Custom_mergedexcess.fits", "Results_MSH_15_52_0_6_1_2_TeV_Custom_mergedexcess.fits", "Results_MSH_15_52_1_2_2_4_TeV_Custom_mergedexcess.fits", "Results_MSH_15_52_gt2_4_TeV_Custom_mergedexcess.fits", "Results_MSH_15_52_std_2d0_ModelPlus_Std_mergedexcess.fits"]



filename_psf      = ["psf_MSH_15_52_0_3_0_6_TeV.fits", "psf_MSH_15_52_0_6_1_2_TeV.fits", "psf_MSH_15_52_1_2_2_4_TeV.fits", "psf_MSH_15_52_gt2_4_TeV.fits", "psf_std.fits"]

# This exposure maps where first generated by the gammapy script then pixel-corrected
filename_expo     = [""] * (nb_pa_analyses)
for a in range(0,nb_pa_analyses):
    filename_expo[a]     = 'exposure_map_' + tags[a] + '_PApix' + '.fits.gz'
    #filename_expo[a] = 'exposure_map_' + tags[a] + '.fits'

# For the X-ray templates of MSH 15-52
filename_thresh   = ["MSH1552_Chandra_v0_4-7keV_thresh_exp8_100x100.fits.gz"] * (nb_pa_analyses)
filename_x        = (["MSH1552_Chandra_vsmo_4-7keV.fits"] * nb_pa_analyses)

filename_slice    = "slices_several.txt" #"slices_perp.txt"


# Fit initial parameters

exposure_amplitude     = [1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-12]
G_comp_ampl            = [1., 1., 1., 0.5, 1.0]
fwhm_init              = 8.

# Analysis pointing coordinates
ra_pt   = 228.5290
dec_pt  = - 59.1575
res_ana = 0.01

# PSR coordinates
ra_psr   = 228.4818
dec_psr  = -59.1358
glon_psr = 320.3208
glat_psr = -1.1619
X_psr    = 201.921740646
Y_psr    = 201.670122119


#Slices to be fit only in an xr conf
if (sliced_selection.find('y') != -1):
    x_slice,y_slice,L_slice,l_slice,dev_slice     =  np.loadtxt(save_path_sliced[analysis] + filename_slice,dtype=float,usecols=(0,1,2,3,4),unpack=True,skiprows=1)

    nb_of_slice   = int(sys.argv[5])

    boxes         = [x_slice[nb_of_slice], y_slice[nb_of_slice], L_slice[nb_of_slice], l_slice[nb_of_slice], dev_slice[nb_of_slice]]




#MultiEnergyBin
path_multiEbins         = "/home/tsirou/Documents/Analyses/MSH_15-52/energy-dependent/MultiEnergyBins/"

filename_Bsteps         = "big_steps/data/energy_excess_b.dat"
filename_Ssteps         = "small_steps/data/energy_excess_s.dat"
filename_Psteps         = "previous/energy_excess_p.dat"
filename_Nsteps         = "new/data/energy_excess_n.dat"
filename_Hsteps         = "h_res/data/energy_excess_h.dat"

path_multiEbins_results = path_multiEbins + "plots/"
