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
# Various variables
#==============================================================================


random_convert_shift =  0  # 1 or 0. Depends on the matrix conventions
display_ds9          = 'y' # yes (y) or no (n) for displaying with ds9 the fit steps
circular_selection   = 'n' # yes (y) or no (n) for selecting a circular region for fitting
bkg_thaw             = 'n' # yes (y) or no (n) for thawing the background amplitude
bkg_amplitude        = 1.0 # usually set to 1.0 but could also be +/- 2%
confidence_levels    = 'y' # yes (y) or no (n) for computing the confidence intervals (time-consuming!)


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
analysis          = 1
# Tag for the different analyses :
#     0    : PA 2.5deg std
#     1    : PA 2.5deg 300 - 600 GeV
#     2    : PA 2.5deg 600 - 900 GeV
#     3    : PA 2.5deg 0.9 - 3.0 TeV
#     4    : PA 2.5deg > 3.0 TeV
###############################################

nb_pa_analyses    = 5

pulsar_then_rebin = 'n' # Check in case one should multiply the x-ray templates before rebinning them --Not fully implemented
gif_norm          = 'y' # For the normalisation of the colorbar in the gifs generation

psf_fudge         = -1  # 0.02 or -1 value in sigma degrees! Has to be converted into pixels. If there is no fudge factor then insert -1

method            = "simplex" # "moncar"

dim_for_a         = '400x400'
# For the hdu FITS extensions for the gamma-like maps
hdu_ext = 1  # 0 for gammapy maps, 1 for PA

# psf_fudge and bkg_amplitude arguments for the systematics estimation
if(analysis == 0 ):
    alpha_factor_pa   = [0.0, 0.0, 1.28, 1.28, 1.27, 1.06, 0.99]
if(analysis == 1 ):
    alpha_factor_pa   = [0.0, 0.0, 1.48, 1.38, 1.37, 1.23, 1.23]
if(analysis == 2 ):
    alpha_factor_pa   = [0.0, 0.0, 1.19, 1.05, 1.01, 0.99, 1.03]
if(analysis == 3 ):
    alpha_factor_pa   = [0.0, 0.0, 1.37, 1.7, 1.9, 1.32, 1.7]
if(analysis == 4 ):
    alpha_factor_pa   = [0.0, 0.0, 0.69, 0.69, 0.57, 0.2, 0.57]

# Display GIF normalisation for the colorbar
cb_extrema              = np.ndarray((nb_pa_analyses, 2))
cb_extrema[0,:]         = [1.0, 4.0]
cb_extrema[1,:]         = [1.0, 2.0]
cb_extrema[2,:]         = [1.0, 2.0]
cb_extrema[3,:]         = [1.0, 3.0]
cb_extrema[4,:]         = [1.0, 2.0]

cb_resids               = np.ndarray((nb_pa_analyses, 2))
cb_resids[0,:]         = [-10.0, 10.0]
cb_resids[1,:]         = [-4.0, 4.0]
cb_resids[2,:]         = [-4.0, 4.0]
cb_resids[3,:]         = [-4.0, 4.0]
cb_resids[4,:]         = [-4.0, 4.0]


#For the exposure_gamma.py script
analysis_energy_low        = [0.3, 0.3, 0.6, 0.9, 3.0]
analysis_energy_high       = [100., 0.6, 0.9, 3.0, 100.]


source            = ["MSH15-52---std","MSH15-52---300_600GeV", "MSH15-52---600_900GeV","MSH15-52---0.9-3.0TeV","MSH15-52---gt3TeV"]
tags              = ["std", "0-3_0-6TeV", "0-6_0-9TeV", "0-9_3-0TeV", "gt3-0TeV"]
fit_component     = ["xr0","xr0-2G","xr","xr-2G","xr-G1","xr-sh","xr-dc"]
#fit_component     = ["xr_psf","xr_psf","xr-2G_psf","xr-G1_psf","xr-sh_psf","xr-dc_psf"]
#fit_component     = ["xr-2G_bkg_d"]



directory         = ["energy-dependent/"] * nb_pa_analyses

path_data         = "/home/tsirou/Documents/Analyses/DATA/"
path_g            = ["/home/tsirou/Documents/Analyses/MSH_15-52/"] * (nb_pa_analyses)
path              = [""] * (nb_pa_analyses)
for n in range(0,nb_pa_analyses):
    path[n]   = path_g[n] + directory[n]

path_x            = "/home/tsirou/Documents/Analyses/Xray_MSH_15-52/"


path_psf          = [p + "psf/" for p in path]
path_template     = [p + "template/" for p in path]
path_expo         = [p + "exposure/" for p in path]
path_runs         = [p + "runs/" for p in path]

save_path         = [p + "Fit/" for p in path]
save_path_gif     = [s + 'GIF/' for s in save_path]
# Path for all the results txt files (global)
results_path      = "/home/tsirou/Documents/Analyses/Results/Energy-dependent/"

filename_runs     = ["Results_MSH_15_52_ModelPlus_Std_mergedrun_list.txt", "Results_MSH_15_52_300_600GeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_600_900GeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_0_9_3TeV_Custom_mergedrun_list.txt", "Results_MSH_15_52_gt3TeV_Custom_mergedrun_list.txt"]

filename_bkg      = ["Results_MSH_15_52_ModelPlus_Std_mergedbackground.fits", "Results_MSH_15_52_300_600GeV_Custom_mergedbackground.fits", "Results_MSH_15_52_600_900GeV_Custom_mergedbackground.fits", "Results_MSH_15_52_0_9_3TeV_Custom_mergedbackground.fits", "Results_MSH_15_52_gt3TeV_Custom_mergedbackground.fits"]
filename_gamma    = ["Results_MSH_15_52_ModelPlus_Std_mergedgamma-like.fits", "Results_MSH_15_52_300_600GeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_600_900GeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_0_9_3TeV_Custom_mergedgamma-like.fits", "Results_MSH_15_52_gt3TeV_Custom_mergedgamma-like.fits"]
filename_excess   = ["Results_MSH_15_52_ModelPlus_Std_mergedexcess.fits", "Results_MSH_15_52_300_600GeV_Custom_mergedexcess.fits", "Results_MSH_15_52_600_900GeV_Custom_mergedexcess.fits", "Results_MSH_15_52_0_9_3TeV_Custom_mergedexcess.fits", "Results_MSH_15_52_gt3TeV_Custom_mergedexcess.fits"]

filename_psf      = ["psf_std.fits", "psf_300_600GeV.fits", "psf_600_900GeV.fits", "psf_0_9_3TeV.fits", "psf_gt3TeV.fits"]

# This exposure maps where first generated by the gammapy script then pixel-corrected by Yves
filename_expo     = [""] * (nb_pa_analyses)
for a in range(0,nb_pa_analyses):
    filename_expo[a]     = 'exposure_map_' + tags[a] + '_PApix' + '.fits.gz'

# For the X-ray templates of MSH 15-52
filename_thresh   = ["MSH1552_Chandra_v0_4-7keV_thresh_exp8_100x100.fits.gz"] * (nb_pa_analyses)
filename_x        = (["MSH1552_Chandra_vsmo_4-7keV.fits"] * nb_pa_analyses)

filename_slice    = ""


# Fit initial parameters

exposure_amplitude     = [1.0e-12, 1.0e-9, 1.0e-10, 1.0e-11,1.0e-12]
G_comp_ampl            = [10., 1., 1., 1., 1.]