import sys
import os
import names

from fitting import fit_1d_psf,fit_2d_psf,fit_bkg,fit_ls5039,fit_msh1552,fit_msh1552_slices,fit_msh1552_xfact_comp
from fitting import rebinning_exp, rebinning_psf,fit_fudged_psf
from results import display_results,quick_results,quick_errors
from astropy import wcs

from xray import xray_map_resize,xray_header
from pulsar import pulsar_mag
from templates import leptonic_bkg,resolution_maps

from display import fits_gif,plot_display,fits_png
from conversion import progress

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy import interpolate


name_of_fits = []

cash         = []
psr_alpha    = []
fits_files   = []
x_ampl       = []
g1           = []


analysis          = names.analysis

hdu               = fits.open(names.path_template[analysis] + names.filename_gamma[analysis])
dim_im            = len(hdu[names.hdu_ext].data)

hdulist           = fits.open(names.path_psf[analysis] + names.filename_psf[analysis])
dim_original_psf  = len(hdulist[names.hdu_ext].data)

hduexp            = fits.open(names.path_expo[analysis] + names.filename_expo[analysis])

factor_psf  = 10
dim_psf     = dim_original_psf / factor_psf

# # Generating the leptonic smoothed background
leptonic_bkg(7.0)

# Background model
bkg = fit_bkg(dim_im,-1)
  

# # Rebinning the PSF image if necessary
rebinning_psf(factor_psf,hdulist,hdu)


if(names.psf_fudge != -1):
    fit_fudged_psf(dim_im)

component_list_index_start   = 0 # len(names.fit_component) - 2
component_list_index_end     = 4 # len(names.fit_component)

for t in range(component_list_index_start, component_list_index_end):


    plt.clf()
    plt.cla()

    cash       = []
    psr_alpha  = []
    fits_files = []
    x_ampl     = []
    g1         = []

    tname = names.fit_component[t] + names.tags[names.analysis]

    # Big steps when varying the alpha parameter
    # steps = 30 #1
    # for f in range(0,steps + 1):
    #     print "\n"
    #     progress(f,steps + 1,tname+'\n')
    #     psr_factor = f * 0.1 #names.alpha_factor_pa[t]
    #
    #     #Rebinning the xray map
    #     xray_map_resize(tname)
    #     #Multiplying it by the distance to the pulsar to a power of psr_factor
    #     pulsar_mag(psr_factor)
    #     #Saving the corresponding header properties in the rebinned X-ray image
    #     xray_header()
    #     # Morphological fit for MSH 15-52 with the X-ray templates
    #     psf,gauss,name_file,component = fit_msh1552_xfact_comp(dim_im,bkg,tname,str(psr_factor),-1)
    #     g1.append(component)
    #     # Loading the CASH stat value for each test
    #     nf,cash_stats    =  np.loadtxt(names.save_path[names.analysis]+"fit_info_"+gauss[:3]+".out",dtype=float,usecols=(0,1),unpack=True,skiprows=1)
    #     cash.append(cash_stats[-1])
    #     psr_alpha.append(psr_factor)
    #     print "\n alpha = ",psr_factor,"CASH = ",cash_stats[-1]
    #     fits_files.append(name_file)
    #
    #     quick_results('no',gauss,name_file,component,psr_factor,cash_stats[-1],-1)
    #
    #
    # #fits_gif(names.save_path,names.save_path_gif,fits_files,cash,psr_alpha,tname,'all_y',g1)
    # fits_gif(names.save_path[analysis],names.save_path_gif[analysis],fits_files,cash,psr_alpha,tname + '_wide','all_n',g1)
    #
    #
    #
    # plt.clf()
    #
    #
    # fig1 = plot_display(20.000,20.000,1000)
    #
    # plt.title(' CASH test for different R$^{a}$',fontsize=20,y=1.02)
    # plt.xlabel('a ',fontsize=22,y=1.02)
    # plt.ylabel('Statistic',fontsize=22,x=1.02)
    #
    # interp = interpolate.interp1d(psr_alpha, cash, kind='cubic')
    # plt.plot(psr_alpha,cash,'b+',markersize=15)
    # plt.plot(psr_alpha,interp(psr_alpha),color='darkgreen')
    # plt.axhline(y=min(cash),color='darkred')
    #
    # plt.axhline(y=min(cash)+1.0,color='orange')
    #
    # for j in range(0,len(cash)):
    #     if( abs(cash[j] - min(cash)) < 0.00005  ):
    #         plt.axvline(x=psr_alpha[j],color='tomato')
    #         plt.annotate("\n\n\n\n\n"+str(psr_alpha[j]),(psr_alpha[j],cash[j]-0.5),fontsize=20,color='SeaGreen')
    #     if( abs(cash[j] - (min(cash) + 1.0)) < 0.005  ):
    #         plt.axvline(x=psr_alpha[j],color='sandybrown')
    #         #plt.annotate("\n\n"+str(psr_alpha[j]),(psr_alpha[j],cash[j]),fontsize=15,color='DarkSeaGreen')
    #
    # plt.savefig(names.save_path[analysis] + 'cash_alpha' + tname + '_wide'+'.pdf',dpi=100)



    plt.clf()
    plt.cla()
    plt.figure()

    cash       = []
    psr_alpha  = []
    fits_files = []
    x_ampl     = []
    g1         = []


    # Small steps when varying the alpha parameter
    # Zoom in for the alpha parameter determination stalled for time related issues


    # steps = 80 #80
    # for f in range(0,steps + 1):
    #     print "\n"
    #     progress(f,steps + 1,tname+'\n')
    #     psr_factor = f * 0.01 + 0.1
    #
    #     xray_map_resize(tname)
    #     #Multiplying it by the distance to the pulsar to a power of psr_factor
    #     pulsar_mag(psr_factor)
    #     #Saving the corresponding header properties in the rebinned X-ray image
    #     xray_header()
    #     # Morphological fit for MSH 15-52 with the X-ray templates
    #     psf,gauss,name_file,component = fit_msh1552_xfact_comp(dim_im,bkg,tname,str(psr_factor),-1)
    #     g1.append(component)
    #     # Loading the CASH stat value for each test
    #     nf,cash_stats    =  np.loadtxt(names.save_path[names.analysis]+"fit_info_"+gauss[:3]+".out",dtype=float,usecols=(0,1),unpack=True,skiprows=1)
    #     cash.append(cash_stats[-1])
    #     psr_alpha.append(psr_factor)
    #     print "\n alpha = ",psr_factor,"CASH = ",cash_stats[-1]
    #     fits_files.append(name_file)
    #
    # fits_gif(names.save_path[analysis],names.save_path_gif[analysis],fits_files,cash,psr_alpha,tname + '_narrow','all_n',g1)
    #
    # plt.clf()
    #
    #
    # fig2 = plot_display(20.000,20.000,1000)
    #
    # plt.title(' CASH test for different R$^{a}$',fontsize=20,y=1.02)
    # plt.xlabel('a ',fontsize=22,y=1.02)
    # plt.ylabel('Statistic',fontsize=22,x=1.02)
    #
    # interp = interpolate.interp1d(psr_alpha, cash, kind='cubic')
    # plt.plot(psr_alpha,cash,'b+',markersize=15)
    # plt.plot(psr_alpha,interp(psr_alpha),color='darkgreen')
    # plt.axhline(y=min(cash),color='darkred')
    #
    # plt.axhline(y=min(cash)+1.0,color='orange')
    #
    # for j in range(0,len(cash)):
    #      if( abs(cash[j] - min(cash)) < 0.00005  ):
    #         plt.axvline(x=psr_alpha[j],color='tomato')
    #         plt.annotate("\n\n\n\n\n"+str(psr_alpha[j]),(psr_alpha[j],cash[j]-0.5),fontsize=20,color='SeaGreen')
    #      if( abs(cash[j] - (min(cash) + 1.0)) < 0.005  ):
    #         plt.axvline(x=psr_alpha[j],color='sandybrown')
    #         #plt.annotate("\n\n"+str(psr_alpha[j]),(psr_alpha[j],cash[j]),fontsize=15,color='DarkSeaGreen')
    #
    # plt.savefig(names.save_path[analysis] + 'cash_alpha_' + tname + '_narrow'+'.pdf',dpi=100)


    plt.clf()
    plt.cla()
    plt.figure()


    # For the results --> plots and confidence intervals

    print "\n"
    progress(t,len(names.fit_component),tname+'\n')
    psr_factor = names.alpha_factor_pa[t]

    #Rebinning the xray map
    xray_map_resize(tname)
    #Multiplying it by the distance to the pulsar to a power of psr_factor
    pulsar_mag(psr_factor)
    #Saving the corresponding header properties in the rebinned X-ray image
    xray_header()
    # Morphological fit for MSH 15-52 with the X-ray templates
    psf,gauss,name_file,component = fit_msh1552_xfact_comp(dim_im,bkg,tname,str(psr_factor),-1)

    # Loading the CASH stat value for each test
    nf,cash_stats    =  np.loadtxt(names.save_path[names.analysis]+"fit_info_"+gauss[:3]+".out",dtype=float,usecols=(0,1),unpack=True,skiprows=1)

    print "\n alpha = ",psr_factor,"CASH = ",cash_stats[-1]

    quick_results('yes',gauss,name_file,component,psr_factor,cash_stats[-1],-1)
    quick_errors('yes',tname, gauss,component, t)

    fits_png(names.save_path[analysis],names.results_path,name_file,cash_stats[-1],psr_factor, tname + '_PA','all_n',component,'n')
    fits_png(names.save_path[analysis], names.results_path, name_file, cash_stats[-1], psr_factor,tname + '_PA', 'all_n', component, 'y')
    fits_png(names.save_path[analysis], names.results_path, name_file, cash_stats[-1], psr_factor, tname + '_PA', 'all_y',
         component, 'n')
    fits_png(names.save_path[analysis], names.results_path, name_file, cash_stats[-1], psr_factor, tname + '_PA', 'all_y',
         component, 'y')

#==============================================================================
# x_slice,y_slice,L_slice,l_slice,dev_slice     =  np.loadtxt(names.path + names.filename_slice,dtype=float,usecols=(0,1,2,3,4),unpack=True,skiprows=1)
# 
# 
# tname = 'negnorm'
# 
# for b in range(0,len(x_slice)):
#     
#     cash       = []
#     psr_alpha  = []
#     fits_files = []
#     x_ampl     = []
#     g1         = []
#     
#     box_slice = [x_slice[b],y_slice[b],L_slice[b],l_slice[b],dev_slice[b]]
#     tname = str(b+1) + tname               
#     steps = 30
# 
#     for f in range(0,steps + 1):
#         print "\n"
#         progress(f,steps + 1,tname+'\n')
#         psr_factor = f * 0.1
#     
#         #Rebinning the xray map
#         xray_map_resize(tname)
#         #Multiplying it by the distance to the pulsar to a power of psr_factor
#         pulsar_mag(psr_factor)    
#         # Morphological fit for MSH 15-52 with the X-ray templates
#         psf,gauss,name_file,component = fit_msh1552_slices(dim_im,bkg,tname,str(psr_factor),box_slice,str(b))
# 
#         g1.append(component)
#         # Loading the CASH stat value for each test
#         nf,cash_stats    =  np.loadtxt(names.save_path+"fit_info_slice_"+str(b)+"_"+gauss[:2]+".out",dtype=float,usecols=(0,1),unpack=True,skiprows=1)
#         cash.append(cash_stats[-1]) 
#         psr_alpha.append(psr_factor)
#         print "\n alpha = ",psr_factor,"CASH = ",cash_stats[-1]
#         fits_files.append(name_file)
# 
#     fits_gif(names.save_path,names.save_path_gif,fits_files,cash,psr_alpha,tname,'all_n',g1)
#  
#     plt.clf()
# 
#     fig1 = plot_display(20.000,20.000,1000) 
# 
#     plt.title('Evolution of the CASH statistic for magnetic component model R$^{a}$',fontsize=20,y=1.02)
#     plt.xlabel('a ',fontsize=20,y=1.02)
#     plt.ylabel('Statistic',fontsize=20,x=1.02)
#  
#     #interp = interpolate.interp1d(psr_alpha, cash, kind='cubic')
#     plt.plot(psr_alpha,cash,'b+',markersize=15)
#     #plt.plot(psr_alpha,interp(psr_alpha),color='darkgreen')
#     plt.axhline(y=min(cash),color='darkred')
# 
# 
#     for j in range(0,len(cash)):
#         if( abs(cash[j] - min(cash)) < 0.0001  ):
#             plt.axvline(x=psr_alpha[j],color='tomato')
#             plt.annotate("\n\n"+str(psr_alpha[j]),(psr_alpha[j],cash[j]),fontsize=15,color='SeaGreen')
# 
#     plt.savefig(names.save_path + 'cash_alpha_neg_norm_slice_'+str(b)+'.pdf',dpi=100)    
#     plt.clf()
#==============================================================================



#
#
# for co in range(0,len(names.fit_component)):
#
#     component  = []
#     cash       = 0.
#     dof        = 0
#
#     tname      = names.fit_component[co] + names.tags[names.analysis]
#     psr_factor = names.alpha_factor_pa[co]
#
#     #Rebinning the xray map
#     xray_map_resize(tname)
#
#     #Multiplying it by the distance to the pulsar to a power of psr_factor
#     pulsar_mag(psr_factor)
#
#     #Saving the corresponding header properties in the rebinned X-ray image
#     xray_header()
#
#     psf,gauss,name_file,component = fit_msh1552_xfact_comp(dim_im,bkg,tname,str(psr_factor),-1)
#
#     nf, cash_stats = np.loadtxt(names.save_path[names.analysis] + "fit_info_" + gauss[:3] + ".out", dtype=float,usecols=(0, 1), unpack=True, skiprows=1)
#     cash  = cash_stats[-1]
#
#     with file(names.save_path[names.analysis] + "fit_info_" + gauss[:3] + ".out") as f:
#         one_row = f.readline()
#     dof = len(one_row.split()) - 3
#     f.close()
#
#     #display_results('yesX',psf,gauss,bkg) # not working well ..
#     quick_results('yes',gauss,name_file,component,psr_factor,cash,dof)
#
#     quick_errors('yes',tname,gauss,component)
#
#     print "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
#
# #==============================================================================
#


hdu.close()
#hdulist.close()
hduexp.close()

sys.stdout.write("\n *** DONE! ***\n")