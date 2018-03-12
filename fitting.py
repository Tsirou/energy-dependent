import sys,os
import options
import names


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sherpa.astro import ui as sherpa
#from sherpa.utils import neville,nearest_interp,neville2d,linear_interp


from scipy.optimize import curve_fit
from scipy import interpolate

from conversion import progress
from conversion import convert_pix2wcs
from conversion import convert_deg2timearcs
from conversion import rebin_med
from conversion import rebin_red
from conversion import crop
from conversion import fwhm_sigma
from conversion import radius_percent
from conversion import gaussian
from conversion import gaussian_gaussian
from conversion import triple_gaussian
from conversion import affine

from display import plot_display

from astropy.io import fits


def fit_bkg(dim,i):

    #==============================================================================


    bkg_option     = "gam"
    save_option    = "y"

    #==============================================================================
    #==============================================================================
    #Loading the background image
    sys.stdout.write("\n*** Background ***\n")
    
    # Type of statistics 
    sherpa.set_method(names.method)
    sherpa.set_stat("cash") #set_stat("cstat")

    if( bkg_option.find("tmp") != -1):        

        sherpa.load_image(names.path_template[names.analysis] + names.filename_bkg[names.analysis])

        # Selecting the type of coordinates
        sherpa.set_coord("logical")

        # Fitting with a gaussian
        sherpa.set_source(sherpa.gauss2d.g0)

        sherpa.freeze(g0.xpos)
        sherpa.freeze(g0.ypos)

        g0.xpos   =   dim/2. + 0.5
        g0.ypos   =   dim/2. + 0.5
        g0.fwhm   =   dim/2.
    
        sherpa.fit()
        #sherpa.image_fit()

        sherpa.thaw(g0.xpos)
        sherpa.thaw(g0.ypos)
    
        sherpa.fit()
        #sherpa.image_fit()
    
        sherpa.freeze(g0.xpos)
        sherpa.freeze(g0.ypos)

        if(save_option.find("y") != -1):
            sherpa.save_model(names.save_path[names.analysis] + names.filename_bkg[names.analysis][:-5]+"_model.fits",clobber=True)
            sherpa.save_resid(names.save_path[names.analysis] + names.filename_bkg[names.analysis][:-5]+"_resid.fits",clobber=True)
        
    if( bkg_option.find("gam") != -1): 
        if(i < 0):
            sherpa.load_table_model("bkg_mod",names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_smoothed_bkg_' + '7.0' + 'sigma.fits')
        else :
            sherpa.load_table_model("bkg_mod",names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_smoothed_bkg_' + str(i) + 'sigma.fits')
        sherpa.show_source()
        
    return bkg_option    
#==============================================================================


def fit_1d_psf():
    
    hdulist                        = fits.open(names.path_psf + names.filename_psf)
    hdulist.info()

    factor_pix2deg                 = abs(hdulist[1].header['cdelt1'])

    all_sum,sum_68,radius_it = radius_percent(68.3,hdulist[1].data,factor_pix2deg)

    print "\nDimension of the PSF image :",hdulist[1].data.shape

    print "*** Theta68 = ",radius_it,"arcmin with ",sum_68," for ",all_sum * 68.3 / 100.," out of ",all_sum  

    dim = len(hdulist[1].data) / 2
    x   = np.arange(0,len(hdulist[1].data))
    y   = np.arange(-int(radius_it / factor_pix2deg / 60.),int(radius_it / factor_pix2deg / 60.))
            
    fig = plt.figure(figsize=(20, 20)) 
    gs  = gridspec.GridSpec(2, 1, height_ratios=[6, 2]) 
    ax0 = plt.subplot(gs[0])

    ax0.set_title('$\Theta$ vs radius for the PSF image',fontsize=20,y=1.02)
    ax0.set_xlabel('$\Theta$ (in arcmin)',fontsize=20,y=1.02)
    ax0.set_ylabel('PSF amplitude',fontsize=20,x=1.02)

    ax0.plot((x-float(dim))*factor_pix2deg *60,hdulist[1].data[dim][x],'g+',markersize=8,label="H.E.S.S. PSF",alpha=0.7)
    ax0.fill_between((y)*factor_pix2deg *60,0,hdulist[1].data[dim][dim+y],color='limegreen',alpha=0.2,label="68.3% inside "+str(radius_it)+" arcmin")

    # Triple Gaussian Fit
    popt,pcov = curve_fit(triple_gaussian,(x-float(dim))*factor_pix2deg *60,hdulist[1].data[dim][x],maxfev=3000)
    ax0.plot((x-dim)*factor_pix2deg *60, triple_gaussian((x-dim)*factor_pix2deg *60, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]), color="darkred",label="3 Gaussians fit",linewidth=1)            
    ax0.plot((x-dim)*factor_pix2deg *60, gaussian((x-dim)*factor_pix2deg *60, popt[0],popt[1],popt[4]), color="gold",linewidth=1)            
    ax0.plot((x-dim)*factor_pix2deg *60, gaussian((x-dim)*factor_pix2deg *60, popt[0],popt[2],popt[5]), color="orange",alpha=0.5,linewidth=1)            
    ax0.plot((x-dim)*factor_pix2deg *60, gaussian((x-dim)*factor_pix2deg *60, popt[0],popt[3],popt[6]), color="orange",linewidth=1)

    ax0.legend(fontsize=15,numpoints=1,loc=1)

    # Double Gaussian Fit
    popt2,pcov2 = curve_fit(gaussian_gaussian,(x-dim)*factor_pix2deg *60,hdulist[1].data[dim][x])
    
    # Simple Gaussian Fit
    popt1,pcov1 = curve_fit(gaussian,(x-dim)*factor_pix2deg *60,hdulist[1].data[dim][x])
    
    # Residuals
    ax1 = plt.subplot(gs[1])        
    
    ax1.plot((x-dim)*factor_pix2deg *60, (hdulist[1].data[dim][x])-triple_gaussian((x-dim)*factor_pix2deg *60, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]), \
         color='red',label="Residuals for 3G fit",linewidth=2)
    ax1.plot((x-dim)*factor_pix2deg *60, (hdulist[1].data[dim][x])-gaussian_gaussian((x-dim)*factor_pix2deg *60, popt2[0],popt2[1],popt2[2],popt2[3],popt2[4]), \
         color='royalblue',label="Residuals for 2G fit",linewidth=2)
    ax1.plot((x-dim)*factor_pix2deg *60, (hdulist[1].data[dim][x])-gaussian((x-dim)*factor_pix2deg *60, popt1[0],popt1[1],popt1[2]), \
         color='purple',label="Residuals for 1G fit",linewidth=2)

    ax1.axhline(color="k", zorder=-1)
    plt.legend(fontsize=15,numpoints=1,loc=1)  
        
    plt.savefig(names.path_psf+'theta_psf_3G'+'.pdf',dpi=100)
    plt.clf()    
    
    hdulist.close()
    
    return 
    
def rebinning_psf(factor,hdulist,hdu):


    if ( str(hdulist[names.hdu_ext].data.shape) != (str(hdu[names.hdu_ext].data.shape)) ):

        im                           = hdulist[names.hdu_ext].data
        print "before rebin or crop : ",len(im)

        im                           = crop(im,490)
        print "before rebin after crop : ",len(im)
        im                           = rebin_med( im, factor )
        print "after rebin and crop : ",len(im)
        
        hdulist[names.hdu_ext].data              = im
        
        hdulist[names.hdu_ext].header['bitpix']  = hdu[names.hdu_ext].header['bitpix']
        hdulist[names.hdu_ext].header['crval1']  = hdu[names.hdu_ext].header['crval1']
        hdulist[names.hdu_ext].header['crval2']  = hdu[names.hdu_ext].header['crval2']
        hdulist[names.hdu_ext].header['crpix1']  = hdu[names.hdu_ext].header['crpix1'] / 8
        hdulist[names.hdu_ext].header['crpix2']  = hdu[names.hdu_ext].header['crpix2'] / 8
        hdulist[names.hdu_ext].header['cdelt1']  = hdu[names.hdu_ext].header['cdelt1']
        hdulist[names.hdu_ext].header['cdelt2']  = hdu[names.hdu_ext].header['cdelt2']
        hdulist[names.hdu_ext].header['ctype1']  = hdu[names.hdu_ext].header['ctype1']
        hdulist[names.hdu_ext].header['ctype2']  = hdu[names.hdu_ext].header['ctype2']


        hdulist.writeto(names.path_psf[names.analysis] + 'rebinned_' + names.filename_psf[names.analysis],clobber=True)
    hdulist.close()
    
    return

def rebinning_exp(factor,hdulist,hdu):

    if ( str(hdulist[0].data.shape) != (str(hdu[1].data.shape)) ):
        im                           = hdulist[0].data

        print "before rebin or crop : ",len(im)
        im                           = rebin_red( im, factor )
        print "after rebin          : ",len(im)
        im                           = crop(im,400)
        print "before rebin and crop : ",len(im)
        
        hdulist[1].data              = im
        hdulist[1].header['crval1']  = hdu[1].header['crval1']
        hdulist[1].header['crval2']  = hdu[1].header['crval2']
        hdulist[1].header['crpix1']  = hdu[1].header['crpix1']
        hdulist[1].header['crpix2']  = hdu[1].header['crpix2']
        hdulist[1].header['cdelt1']  = hdu[1].header['cdelt1']
        hdulist[1].header['cdelt2']  = hdu[1].header['cdelt2']
        hdulist[1].header['ctype1']  = hdu[1].header['ctype1']
        hdulist[1].header['ctype2']  = hdu[1].header['ctype2']

    hdulist.writeto(names.path + 'rebinned_' + names.filename_expo, clobber=True)
    hdulist.close()
    
    return

def fit_2d_psf(dim_psf):
    sherpa.load_image(names.path_psf + "rebinned" + names.filename_psf)
    
    # Fitting the PSF in order to have some info about it
    sherpa.set_stat("leastsq")

    sherpa.set_source(sherpa.gauss2d.psf0 + sherpa.gauss2d.psf1 + sherpa.gauss2d.psf2)    

    psf0.xpos   =  dim_psf / 2.
    psf0.ypos   =  dim_psf / 2.

    psf1.xpos   =  dim_psf / 2.
    psf1.ypos   =  dim_psf / 2.
    
    psf2.xpos   =  dim_psf / 2.
    psf2.ypos   =  dim_psf / 2.
        
    sherpa.freeze(psf0.xpos)
    sherpa.freeze(psf0.ypos)

    sherpa.freeze(psf1.xpos)
    sherpa.freeze(psf1.ypos)

    sherpa.freeze(psf2.xpos)
    sherpa.freeze(psf2.ypos)
    
    sherpa.fit()
    #sherpa.image_fit()
    
    sherpa.thaw(psf0.xpos)
    sherpa.thaw(psf0.ypos)
    sherpa.thaw(psf1.xpos)
    sherpa.thaw(psf1.ypos)   
    sherpa.thaw(psf2.xpos)
    sherpa.thaw(psf2.ypos)   

    sherpa.fit()
    #sherpa.image_fit()
    
    sherpa.freeze(psf0.xpos)
    sherpa.freeze(psf0.ypos)
    
    sherpa.fit()
    #sherpa.image_fit()

    sherpa.freeze(psf1.xpos)
    sherpa.freeze(psf1.ypos)

    sherpa.fit()
    #sherpa.image_fit()
    
    sherpa.freeze(psf2.xpos)
    sherpa.freeze(psf2.ypos)

    sys.stdout.write("\n--- PSF characteristics ---\n")
    sherpa.fit(outfile=save_path+"fit_" + names.filename_psf[:-5] + ".out", clobber=True)
    #sherpa.image_fit()
    sys.stdout.write("\n---------------------------\n")

    return

#==============================================================================
# MSH 15-52 with 1, 2 or 3 Gaussians
#==============================================================================
def fit_msh1552(dim,bkg_option,i):

    # Convoluted PSF model
    sherpa.load_psf( "hPSF2", names.path_psf[names.analysis] + "rebinned_" + names.filename_psf[names.analysis])
    sherpa.set_psf(hPSF2)
    print "\n - - - - - - - - - - - - \n"
    sherpa.show_psf()
    print "\n - - - - - - - - - - - - \n"
    print(hPSF2)
    # sherpa.image_psf()
    # sherpa.set_psf(hPSF2)
    # sherpa.show_kernel()

    #==============================================================================
    #==============================================================================
    # Exposure 
    sherpa.load_table_model("exposure",names.path_expo[names.analysis] + names.filename_expo[names.analysis])
    print "Exposure :\n",exposure
    sherpa.show_source()

    #==============================================================================
    #==============================================================================
        
    # Loading the gamma-like candidates image
    sys.stdout.write("\n*** Gamma-like candidates ***\n")
    sherpa.load_image(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
    
    sherpa.set_method(names.method)
    sherpa.set_stat("cash")

    # Using the parameter for the background and adding a gaussian fit for the source

    if (i < 0):
        psf        = options.use_psf()
        psf_gauss  = options.nb_gauss()
    else :
        psf        = "y"
        psf_gauss  = "G3"

    if(psf.find("y") != -1 and bkg_option.find("tmp") != -1 ):
        if(psf_gauss.find("G3") != -1):
            sherpa.set_full_model(g0 + hPSF2((sherpa.gauss2d.g1 + sherpa.gauss2d.g2 + sherpa.gauss2d.g3) * exposure))
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(g0 + hPSF2(sherpa.gauss2d.g1 * exposure))
        if(psf_gauss.find("G2") != -1):
            sherpa.set_full_model(g0 + hPSF2((sherpa.gauss2d.g1 + sherpa.gauss2d.g2) *exposure ))

    if(psf.find("n") != -1 and bkg_option.find("tmp") != -1):
        if(psf_gauss.find("G3") != -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1 + sherpa.gauss2d.g2 + sherpa.gauss2d.g3)
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1)
        if(psf_gauss.find("G2") != -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1 + sherpa.gauss2d.g2)

    if(psf.find("y") != -1 and bkg_option.find("gam") != -1 ):
        if(psf_gauss.find("G3") != -1):
            sherpa.set_full_model(bkg_mod + hPSF2( (sherpa.gauss2d.g1 + sherpa.gauss2d.g2 + sherpa.gauss2d.g3) * exposure))
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(bkg_mod + hPSF2(sherpa.gauss2d.g1 * exposure))
        if(psf_gauss.find("G2") != -1):
            sherpa.set_full_model(bkg_mod + hPSF2( (sherpa.gauss2d.g1 + sherpa.gauss2d.g2) * exposure ))
        
    if(psf.find("n") != -1 and bkg_option.find("gam") != -1):
        if(psf_gauss.find("G3") != -1):
            sherpa.set_full_model(bkg_mod + sherpa.gauss2d.g1 + sherpa.gauss2d.g2 + sherpa.gauss2d.g3)
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(bkg_mod + sherpa.gauss2d.g1)
        if(psf_gauss.find("G2") != -1):
            sherpa.set_full_model(bkg_mod + sherpa.gauss2d.g1 + sherpa.gauss2d.g2)              


    #==============================================================================
    #==============================================================================

    exposure.ampl   =   1.0e-12
    sherpa.freeze(exposure.ampl)

    bkg_mod.ampl = 1.0
    sherpa.freeze(bkg_mod.ampl)
        
    # 3G 2D fit
    if(psf_gauss.find("G3") != -1):        
    
        # A gaussian on the left side
        g1.xpos   =   dim/2. - 2.
        g1.ypos   =   dim/2. - 1.5
        sherpa.set_par(g1.ampl, min = 0.)

        # A gaussian on the right side
        g2.xpos   =   dim/2. + 3.5
        g2.ypos   =   dim/2. + 3.0
        sherpa.set_par(g2.ampl, min = 0.)        
        
        # A gaussian at the center
        g3.xpos   =   dim/2.
        g3.ypos   =   dim/2.
        g3.ampl   =   3
        sherpa.set_par(g3.ampl, min = 0.)        
        sherpa.set_par(g3.fwhm, min = 10.)        
        
        
        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)
        sherpa.freeze(g2.xpos)
        sherpa.freeze(g2.ypos)
        sherpa.freeze(g3.xpos)
        sherpa.freeze(g3.ypos)
        sherpa.freeze(g3.ampl)
      
        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g1.xpos)
        sherpa.thaw(g1.ypos)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g2.xpos)
        sherpa.thaw(g2.ypos)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g3.ampl)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g3.xpos)
        sherpa.thaw(g3.ypos)
        
        
        
    # 1G 2D fit        
    if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):

        g1.xpos   =   dim/2. + 0.5
        g1.ypos   =   dim/2. + 0.5
        g1.ellip  =   0
        g1.theta  =   0

        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)
        sherpa.freeze(g1.ellip)
        sherpa.freeze(g1.theta)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g1.xpos)
        sherpa.thaw(g1.ypos)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g1.ellip)
        sherpa.thaw(g1.theta)


    # 2G 2D fit
    if(psf_gauss.find("G2") != -1):

        # A gaussian on the left side
        g1.xpos   =   dim/2. - 2.
        g1.ypos   =   dim/2. - 1.5

        # A gaussian on the right side
        g2.xpos   =   dim/2. + 3.5
        g2.ypos   =   dim/2. + 3.0

        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)

        sherpa.freeze(g2.xpos)
        sherpa.freeze(g2.ypos)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g1.xpos)
        sherpa.thaw(g1.ypos)

        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        sherpa.thaw(g2.xpos)
        sherpa.thaw(g2.ypos)    

    sys.stdout.write("\n###########################################\n")
    if (i == -1 ):
        sherpa.fit(outfile=names.save_path[names.analysis] + "fit_info_"+psf_gauss[:2]+".out", clobber=True)
    else :
        sherpa.fit(outfile=names.save_path[names.analysis] + "fit_"+str(i)+"_"+psf_gauss[:2]+".out", clobber=True)

    print "\n"
    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()
    print "\n"
    sherpa.show_model()
    sys.stdout.write("\n###########################################\n")    
    
    # Saving the model 

    if(psf.find("y") != -1):
        sherpa.save_model(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + "_model.fits",clobber=True)
        sherpa.save_resid(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + "_resid.fits",clobber=True)
        fits_filename              = names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + "_model.fits"

    if(psf.find("n") != -1):
        sherpa.save_model(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + "_nopsf_model.fits",clobber=True)
        sherpa.save_resid(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + "_nopsf_resid.fits",clobber=True)
        fits_filename              = names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + "_nopsf_model.fits"

    if (i == -1):
        sys.stdout.write("\n__________________________________________\nConfidence intervals\nCompute? y/n :\n")
        confidence   = sys.stdin.readline()

        if( confidence.find("y") != -1):
            sherpa.conf()
            temporary = sys.stdout
            sys.stdout = open(names.save_path + "conf_" + psf_gauss[:2] + ".txt", 'w')
            print(sherpa.get_conf_results())
            sys.stdout.close()
            sys.stdout = temporary
        sys.stdout.write("\n__________________________________________\n")
        
    hdu                            = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

    factor_pix2deg_gamma           = abs(hdu[names.hdu_ext].header['cdelt1'])
    all_sum_g,sum_68_g,radius_it_g = radius_percent(68.3,hdu[names.hdu_ext].data,factor_pix2deg_gamma)
    print "Theta68 = ",radius_it_g,"arcmin with ",sum_68_g," for ",all_sum_g * 68.3 / 100.," out of ",all_sum_g  

    hdu.close()

    if(psf_gauss.find("G3") != -1):
        sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             

        sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
        alpha_gauss2,delta_gauss2   = convert_pix2wcs(fits_filename, g2.xpos.val, g2.ypos.val)
        convert_deg2timearcs(alpha_gauss2,delta_gauss2)             

        sys.stdout.write("\n Coordinates of the gaussian at the center :\n")
        alpha_gauss3,delta_gauss3   = convert_pix2wcs(fits_filename, g3.xpos.val, g3.ypos.val)
        convert_deg2timearcs(alpha_gauss3,delta_gauss3)             

        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (left gaussian)  : ",r_intrinsic," arcmin"

        r_intrinsic         = fwhm_sigma(g2.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (right gaussian) : ",r_intrinsic," arcmin"    
    
        r_intrinsic         = fwhm_sigma(g3.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (center)         : ",r_intrinsic," arcmin"
    
    if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
        sys.stdout.write("\n Coordinates of the gaussian :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             

        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation : ",r_intrinsic," arcmin"


    if(psf_gauss.find("G2") != -1):
        sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             
    
        sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
        alpha_gauss2,delta_gauss2   = convert_pix2wcs(fits_filename, g2.xpos.val, g2.ypos.val)
        convert_deg2timearcs(alpha_gauss2,delta_gauss2)             
    
        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (left gaussian)  : ",r_intrinsic," arcmin"

        r_intrinsic         = fwhm_sigma(g2.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (right gaussian) : ",r_intrinsic," arcmin"    
    
    return psf,psf_gauss 
    
def fit_ls5039(dim):

    # Convoluted PSF model
    sherpa.load_psf( "hPSF2", names.path_psf + "rebinned_" + names.filename_psf)
    sherpa.set_psf(hPSF2)
    print "\n - - - - - - - - - - - - \n"
    sherpa.show_psf()
    print "\n - - - - - - - - - - - - \n"
    print(hPSF2)
    sherpa.image_psf()

    sherpa.set_psf(hPSF2) 
    sherpa.show_kernel()

    #==============================================================================
    
    # Loading the gamma-like candidates image
    sys.stdout.write("\n*** Gamma-like candidates ***\n")
    sherpa.load_image(names.path + names.filename_gamma)

    sherpa.set_stat("cash")

    # Using the parameter for the background and adding a gaussian fit for the source

    psf   = options.use_psf()
    if(psf.find("y") != -1):
        psf_gauss = "G1"
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(g0 + hPSF2(sherpa.gauss2d.g1))
        
    if(psf.find("n") != -1):
        psf_gauss = options.nb_gauss()        
        if(psf_gauss.find("G3") != -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1 + sherpa.gauss2d.g2 + sherpa.gauss2d.g3)
        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1)
        if(psf_gauss.find("G2") != -1):
            sherpa.set_full_model(g0 + sherpa.gauss2d.g1 + sherpa.gauss2d.g2)      

    if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
        # A gaussian
        g1.xpos   =   dim/2.
        g1.ypos   =   dim/2.
        g1.fwhm   =   1.84

        sherpa.freeze(g1.fwhm)
        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)
    
        sherpa.fit()
        #sherpa.image_fit()
        
        sherpa.thaw(g1.xpos)
        sherpa.thaw(g1.ypos)

        sherpa.fit()
        #sherpa.image_fit()

        likelihood_profile             = options.cash_profile()
        
        if(likelihood_profile.find("n") != -1):     
            sherpa.thaw(g1.fwhm)

        if(likelihood_profile.find("y") != -1):             
            fwhm_val  = []
            sigma_val = []
            cash_val  = []

            sherpa.freeze(g1.xpos)
            sherpa.freeze(g1.ypos)

            for i in range(0,450):
                progress(i,450,'')
            
                nfev           = []
                statistic_cash = []

                g1.fwhm  = 0.6 + i * 0.01
                fwhm_val.append(0.6 + i * 0.01)   

                #==============================================================================
                # g1.fwhm  = 0.2 + i * 0.01
                # fwhm_val.append(0.2 + i * 0.01)    
                #==============================================================================
                sherpa.freeze(g1.fwhm)

                sys.stdout.write("\n###########################################\n")
                sherpa.fit(outfile=names.save_path + "fit_info_ls5039_"+str(i)+".out", clobber=True)
                print "\n"
                #sherpa.image_fit()
                print "\n"
                sherpa.show_model()
                sys.stdout.write("\n###########################################\n")    
    
    
                nfev,statistic_cash    =  np.loadtxt(names.save_path+"fit_info_ls5039_"+str(i)+".out",dtype=float,usecols=(0,1),unpack=True,skiprows=1)
                cash_val.append(statistic_cash[-1])

            f  = open(names.save_path + "ls5039_cash_vs_fwhm.txt",'w')    
            f.write("fwhm CASH\n")
            for t in range(0,len(cash_val)):
                f.write(str(fwhm_val[t]) + " " + str(cash_val[t])+"\n")
            f.close()
    
            sys.stdout.write("\n###########################################\nTHAWED FWHM\n")

            sherpa.thaw(g1.fwhm)
            sherpa.freeze(g1.xpos)
            sherpa.freeze(g1.ypos)

        sherpa.fit()
        print "\n"
        #sherpa.image_fit()
        print "\n"
        sherpa.show_model()
        sys.stdout.write("\n###########################################\n") 
 
        hdu                            = fits.open(names.path + names.filename_gamma)
    
        factor_pix2deg_gamma           = abs(hdu[1].header['cdelt1'])
        all_sum_g,sum_68_g,radius_it_g = radius_percent(68.3,hdu[1].data,factor_pix2deg_gamma)
        print "Theta68 = ",radius_it_g,"arcmin with ",sum_68_g," for ",all_sum_g * 68.3 / 100.," out of ",all_sum_g  

        hdu.close()
    
        if(likelihood_profile.find("y") != -1): 
            fig = plot_display(20.000, 20.000 , 1000)
            
            plt.clf()    

            plt.title('Evolution of the CASH statistic for given FWHM values',fontsize=20,y=1.02)
            plt.xlabel('$\sigma_{intrinsic}$ (arcsec)',fontsize=20,y=1.02)
            plt.ylabel('Statistic',fontsize=20,x=1.02)

            for i in range(0,len(fwhm_val)):
                sigma_val.append(fwhm_sigma(fwhm_val[i], factor_pix2deg_gamma)*60.)

            plt.plot(sigma_val, cash_val, color="darkblue", marker="+",label="Cash statistic",markersize=3)            

            #==============================================================================
            # confidence_level = [0.,1.0,2.3,3.53,4.22]
            # colours          = ["darkgreen","purple","magenta","red","orange"]
            # epsilon          = [1.0e-30,0.001,0.1,0.1,0.1]
            #==============================================================================
            confidence_level = [0.,1.0,4.0,9.0]
            colours          = ["darkgreen","purple","red","orange"]
            epsilon          = [1.0e-30,0.001,0.2,0.1 ]
            
            for i in range(0,len(confidence_level)):
                inter            = []
                #plt.axhline(y=min(cash_val)+confidence_level[i],color=colours[i],label=str(i)+" dof")
                plt.axhline(y=min(cash_val)+confidence_level[i],color=colours[i],label=str(i)+" $\sigma$")
    
                for j in range(0,len(cash_val)):
                    if(abs(cash_val[j] - (min(cash_val)+confidence_level[i])) < epsilon[i]  ):
                        inter.append(j)
                        print "For ",i," dof : ",cash_val[j],fwhm_val[j]

                plt.axvline(x=sigma_val[min(inter)],color=colours[i])
                plt.axvline(x=sigma_val[max(inter)],color=colours[i])
        
                plt.annotate("  "+str(sigma_val[min(inter)])[:4],(sigma_val[min(inter)],cash_val[min(inter)]),fontsize=20,color=colours[i])
                plt.annotate("  "+str(sigma_val[max(inter)])[:4],(sigma_val[max(inter)],cash_val[max(inter)]),fontsize=20,color=colours[i])

                plt.legend(fontsize=15,numpoints=1,loc=1)          

            plt.legend(fontsize=15,numpoints=1,loc=1)          
            plt.savefig(names.save_path + 'ls5039_cash_arcsec'+'.pdf',dpi=100)
            plt.clf()   

    if(psf_gauss.find("G2") != -1):        
        # 1G gaussian
        g1.xpos   =   dim/2. + 0.5
        g1.ypos   =   dim/2. + 0.5
        sherpa.set_par(g1.ampl,min=0.)

        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)

        # 2G gaussian
        g2.xpos   =   dim/2. + 0.5
        g2.ypos   =   dim/2. + 0.5
        sherpa.set_par(g2.ampl,min=0.)
        
        sherpa.freeze(g2.xpos)
        sherpa.freeze(g2.ypos)
        
        sherpa.fit()
        #sherpa.image_fit()
        
        sherpa.thaw(g1.xpos)
        sherpa.thaw(g1.ypos)

        sherpa.fit()
        #sherpa.image_fit()

        sherpa.thaw(g2.xpos)
        sherpa.thaw(g2.ypos)

        sherpa.fit()
        #sherpa.image_fit()

    if(psf_gauss.find("G3") != -1):        
        # 1G gaussian
        g1.xpos   =   dim/2. + 0.5
        g1.ypos   =   dim/2. + 0.5
        sherpa.set_par(g1.ampl,min=0.)

        sherpa.freeze(g1.xpos)
        sherpa.freeze(g1.ypos)

        # 2G gaussian
        g2.xpos   =   dim/2. + 0.5
        g2.ypos   =   dim/2. + 0.5
        sherpa.set_par(g2.ampl,min=0.)
        
        sherpa.freeze(g2.xpos)
        sherpa.freeze(g2.ypos)

        # 3G gaussian
        g3.xpos   =   dim/2. + 0.5
        g3.ypos   =   dim/2. + 0.5
        sherpa.set_par(g3.ampl,min=0.)

        sherpa.freeze(g3.xpos)
        sherpa.freeze(g3.ypos)
        
        sherpa.fit()
        #sherpa.image_fit()
        
   
    #==============================================================================
    #==============================================================================
    # Saving the model 

    if(psf.find("y") != -1):
        sherpa.save_model(names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_model.fits",clobber=True)
        sherpa.save_resid(names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_resid.fits",clobber=True)
        fits_filename              = names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_model.fits"

    if(psf.find("n") != -1):
        sherpa.save_model(names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_nopsf_model.fits",clobber=True)
        sherpa.save_resid(names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_nopsf_resid.fits",clobber=True)
        fits_filename              = names.save_path + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + "_nopsf_model.fits"

    sys.stdout.write("\n__________________________________________\nConfidence intervals\nCompute? y/n :\n")
    confidence   = sys.stdin.readline()
    if( confidence.find("y") != -1):
        sherpa.conf()
        temporary = sys.stdout
        sys.stdout = open(names.save_path + "conf_ls5039_"+psf_gauss[:2]+".txt", 'w')
        print(sherpa.get_conf_results())
        sys.stdout.close()
        sys.stdout = temporary        
    sys.stdout.write("\n__________________________________________\n")

    hdu                            = fits.open(names.path_template + names.filename_gamma)

    factor_pix2deg_gamma           = abs(hdu[1].header['cdelt1'])
    all_sum_g,sum_68_g,radius_it_g = radius_percent(68.3,hdu[1].data,factor_pix2deg_gamma)
    print "Theta68 = ",radius_it_g,"arcmin with ",sum_68_g," for ",all_sum_g * 68.3 / 100.," out of ",all_sum_g  

    hdu.close()


    if(psf_gauss.find("G3") != -1):
        sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             

        sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
        alpha_gauss2,delta_gauss2   = convert_pix2wcs(fits_filename, g2.xpos.val, g2.ypos.val)
        convert_deg2timearcs(alpha_gauss2,delta_gauss2)             

        sys.stdout.write("\n Coordinates of the gaussian at the center :\n")
        alpha_gauss3,delta_gauss3   = convert_pix2wcs(fits_filename, g3.xpos.val, g3.ypos.val)
        convert_deg2timearcs(alpha_gauss3,delta_gauss3)             

        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (left gaussian)  : ",r_intrinsic," arcmin"

        r_intrinsic         = fwhm_sigma(g2.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (right gaussian) : ",r_intrinsic," arcmin"    
    
        r_intrinsic         = fwhm_sigma(g3.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (center)         : ",r_intrinsic," arcmin"
    
    if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
        sys.stdout.write("\n Coordinates of the gaussian :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             

        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation : ",r_intrinsic," arcmin"


    if(psf_gauss.find("G2") != -1):
        sys.stdout.write("\n Coordinates of the gaussian on the left :\n")
        alpha_gauss1,delta_gauss1   = convert_pix2wcs(fits_filename, g1.xpos.val, g1.ypos.val)
        convert_deg2timearcs(alpha_gauss1,delta_gauss1)             
    
        sys.stdout.write("\n Coordinates of the gaussian on the right :\n")
        alpha_gauss2,delta_gauss2   = convert_pix2wcs(fits_filename, g2.xpos.val, g2.ypos.val)
        convert_deg2timearcs(alpha_gauss2,delta_gauss2)             
    
        sys.stdout.write("\n Source morphology \n")
        r_intrinsic         = fwhm_sigma(g1.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (left gaussian)  : ",r_intrinsic," arcmin"

        r_intrinsic         = fwhm_sigma(g2.fwhm.val,factor_pix2deg_gamma)
        print "\n Intrinsic standart deviation (right gaussian) : ",r_intrinsic," arcmin"    
    
    return psf,psf_gauss


def fit_fudged_psf(dim):

    psf_hdu = fits.open(names.path_psf[names.analysis] + "rebinned_" + names.filename_psf[names.analysis])

    # Input value in sigma can cnvert it into the actual size
    fudge = names.psf_fudge * (2. * np.sqrt(2. * np.log(2.)))

    # Converting degrees to pixels by reading the binsize from the PSF image header
    psf_binsize = abs(psf_hdu[names.hdu_ext].header['cdelt1'])
    psf_hdu.close()

    fudge = fudge / psf_binsize

    print "Fudge of ", names.psf_fudge, " sigma (deg) --> ", fudge, " pixels"

    sherpa.load_psf("hPSF2", names.path_psf[names.analysis] + "rebinned_" + names.filename_psf[names.analysis])
    sherpa.set_psf(hPSF2)

    sherpa.load_image(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

    sherpa.set_full_model(hPSF2(sherpa.gauss2d.gpsf))

    gpsf.xpos = dim / 2. + 0.5
    gpsf.ypos = dim / 2. + 0.5

    gpsf.ellip = 0.
    gpsf.theta = 0.

    gpsf.fwhm = fudge

    sherpa.freeze(gpsf.xpos)
    sherpa.freeze(gpsf.ypos)
    sherpa.freeze(gpsf.ellip)
    sherpa.freeze(gpsf.theta)
    sherpa.freeze(gpsf.fwhm)

    gpsf.ampl = 1.0

    sherpa.fit()

    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()

    sherpa.thaw(gpsf.ampl)

    sherpa.fit()

    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()


    sys.stdout.write("\n Saving the fudged psf ... \n")

    sherpa.save_model(1,names.path_psf[names.analysis] + names.filename_psf[names.analysis][:-5] + "_fudge_model.fits",clobber=True)
    sherpa.save_resid(1,names.path_psf[names.analysis] + names.filename_psf[names.analysis][:-5] + "_fudge_resid.fits",clobber=True)

    sherpa.show_model()

    return

#==============================================================================
# MSH 15-52 with the X-ray templates + factor + amplitude             
#==============================================================================
def fit_msh1552_xfact_comp(dim,bkg_option,psf_gauss,psr_alpha,xampl):

    sherpa.set_method(names.method)
    # Convoluted PSF model

    # Using the parameter for the background and adding a gaussian fit for the source

    psf       = "y"
    psf_gauss = "X" + psf_gauss

    if (psf.find("y") != -1 and names.psf_fudge != -1):

        hdu_psf_fudge = fits.open(names.path_psf[names.analysis] + names.filename_psf[names.analysis][:-5] + "_fudge_model.fits")
        im_psf         = hdu_psf_fudge[0].data

        hdu_psf_fudge[0].data = crop(im_psf,49) # Carefull! this is for a 400x400 pix2 analysis

        hdu_psf_fudge.flush()
        hdu_psf_fudge.close()


        sherpa.load_psf("hPSF2",names.path_psf[names.analysis] + names.filename_psf[names.analysis][:-5] + "_fudge_model.fits")

    else:
        sherpa.load_psf("hPSF2", names.path_psf[names.analysis] + "rebinned_" + names.filename_psf[names.analysis])


    # sherpa.set_psf(hPSF2)
    # print "\n - - - - - - - - - - - - \n"
    # sherpa.show_psf()
    # print "\n - - - - - - - - - - - - \n"
    # print(hPSF2)
    # sherpa.image_psf()
    # sherpa.set_psf(hPSF2)
    # sherpa.show_kernel()


    #==============================================================================
    #==============================================================================
    
    # X-ray model    
    sherpa.load_table_model( "xrays_mod", names.path_x + "rebinned_" + names.dim_for_a + ".fits" )
    sherpa.show_source()
    print "Xray model :\n",xrays_mod


    # Exposure
    sherpa.load_table_model("exposure",names.path_expo[names.analysis] + names.filename_expo[names.analysis])
    #sherpa.set_exposure(exposure)
    print "Exposure :\n",exposure
    sherpa.show_source()
    #sherpa.set_data("exposure",exposure)


    #==============================================================================     
    #==============================================================================

    # Loading the gamma-like candidates image
    sys.stdout.write("\n*** Gamma-like candidates ***\n")
    sherpa.load_image(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

    sherpa.set_stat("cash")


    
    if(psf.find("y") != -1):
        if(psf_gauss.find("G1") == -1 and psf_gauss.find("2G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2(xrays_mod * exposure))
        if(psf_gauss.find("G1") != -1 or psf_gauss.find("2G") != -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2((xrays_mod + sherpa.gauss2d.gcomp)) * exposure)
        if(psf_gauss.find("sh") != -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2((xrays_mod + sherpa.shell2d.shell)) * exposure)
        if(psf_gauss.find("dc") != -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2( (xrays_mod + sherpa.disk2d.disc)) * exposure)
        if(psf_gauss.find("Sc") != -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2( (xrays_mod + sherpa.sersic2d.gcomp)) * exposure)



    if(names.circular_selection.find('y') != -1):
        # Selecting a region for the fit procedure
        sherpa.notice2d("circle("+str(dim/2 + 0.5)+","+str(dim/2 + 0.5)+","+str(dim/4)+")")
    if(names.sliced_selection.find('y') != -1):
        sherpa.notice2d("box(" + str(names.boxes[0]) + "," + str(names.boxes[1]) + "," + str(names.boxes[2]) + "," + str(names.boxes[3]) + "," + str(names.boxes[4]) + ")")

                
    if(bkg_option.find("gam") != -1):

        bkg_mod.ampl    =   names.bkg_amplitude
        sherpa.freeze(bkg_mod.ampl)
        
        exposure.ampl   =   names.exposure_amplitude[names.analysis]
        sherpa.freeze(exposure.ampl)

        sherpa.set_par(xrays_mod.ampl, min = 1.)
        xrays_mod.ampl  = 1.0e5 # 1.0e2 #
        
        if(psf_gauss.find("G1") != -1 or psf_gauss.find("2G") != -1):
            gcomp.xpos   =   dim/2. + 0.5
            gcomp.ypos   =   dim/2. + 0.5
            gcomp.ellip  =   0
            gcomp.theta  =   0
            gcomp.ampl   =   names.G_comp_ampl[names.analysis]
            gcomp.fwhm   =   20.
            sherpa.set_par(gcomp.ampl, min=0.)
            sherpa.set_par(gcomp.xpos, min = dim/4., max = 3 * dim/4.)
            sherpa.set_par(gcomp.ypos, min = dim/4., max = 3 * dim/4.)

            sherpa.freeze(gcomp.xpos)
            sherpa.freeze(gcomp.ypos)
            sherpa.freeze(gcomp.ellip)
            sherpa.freeze(gcomp.theta)
            sherpa.freeze(gcomp.ampl)
            sherpa.freeze(gcomp.fwhm)
            

        if(psf_gauss.find("sh") != -1):
            
            shell.xpos   =   dim/2. + 0.5
            shell.ypos   =   dim/2. + 0.5
            shell.ampl   =   names.G_comp_ampl[names.analysis]
            shell.r0     =   10.
            shell.width  =   20.
            sherpa.set_par(shell.ampl, min = 0.)
            sherpa.set_par(shell.xpos, min = 0., max=dim)
            sherpa.set_par(shell.ypos, min = 0., max=dim)
            
            sherpa.freeze(shell.xpos)
            sherpa.freeze(shell.ypos)
            sherpa.freeze(shell.ampl)
            sherpa.freeze(shell.r0)
            sherpa.freeze(shell.width)
            
        if(psf_gauss.find("dc") != -1):
            
            disc.xpos   =   dim/2. + 0.5
            disc.ypos   =   dim/2. + 0.5
            disc.ampl   =   names.G_comp_ampl[names.analysis]
            disc.r0     =   10.

            sherpa.set_par(disc.ampl, min=0.)
            sherpa.set_par(disc.xpos, min = 0., max=dim)
            sherpa.set_par(disc.ypos, min = 0., max=dim)


            sherpa.freeze(disc.xpos)
            sherpa.freeze(disc.ypos)
            sherpa.freeze(disc.ampl)            
            sherpa.freeze(disc.r0)


        if(psf_gauss.find("Sc") != -1):
            gcomp.xpos   =   names.X_psr
            gcomp.ypos   =   names.Y_psr
            gcomp.ellip  =   0
            gcomp.theta  =   0
            gcomp.ampl   =   names.G_comp_ampl[names.analysis]
            gcomp.r0     =   names.fwhm_init
            gcomp.n      =   0.1
            sherpa.set_par(gcomp.ampl, min=0.)
            sherpa.set_par(gcomp.xpos, min = dim/4., max = 3 * dim/4.)
            sherpa.set_par(gcomp.ypos, min = dim/4., max = 3 * dim/4.)

            sherpa.freeze(gcomp.xpos)
            sherpa.freeze(gcomp.ypos)
            sherpa.freeze(gcomp.ellip)
            sherpa.freeze(gcomp.theta)
            sherpa.freeze(gcomp.ampl)
            sherpa.freeze(gcomp.r0)

            sherpa.freeze(gcomp.n)

        sherpa.fit()

        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        if(names.bkg_thaw.find('y') != -1):
            sherpa.thaw(bkg_mod.ampl) # Should it stay as a thawed or a frozen parameter?
            sherpa.fit()
            if(names.display_ds9.find('y') != -1 ):
                sherpa.image_fit()
                                    
        sherpa.thaw(xrays_mod.ampl)   
        
        sherpa.fit()
        if(names.display_ds9.find('y') != -1 ):
            sherpa.image_fit()

        if(xampl > 0.):                            
            xrays_mod.ampl   =   xampl   
            sherpa.freeze(xrays_mod.ampl)

        if(psf_gauss.find("G1") != -1 and psf_gauss.find("2G1") == -1):
        
            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.xpos)
            sherpa.thaw(gcomp.ypos)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.ellip)
            sherpa.thaw(gcomp.theta)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.ampl)

            gcomp.fwhm   =   names.fwhm_init
            sherpa.freeze(gcomp.fwhm)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.fwhm)



        if(psf_gauss.find("2G") != -1):
        
            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.xpos)
            sherpa.thaw(gcomp.ypos)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.ampl)

            gcomp.fwhm   =   8
            sherpa.freeze(gcomp.fwhm)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.fwhm)
            
            
        if(psf_gauss.find("sh") != -1):

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(shell.ampl)
            sherpa.thaw(shell.r0) 


            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(shell.xpos)
            sherpa.thaw(shell.ypos)


            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(shell.width)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.freeze(shell.ampl)            

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(shell.ampl)            

        if(psf_gauss.find("dc") != -1):

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(disc.ampl)
            sherpa.thaw(disc.r0)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(disc.xpos)
            sherpa.thaw(disc.ypos)

        if (psf_gauss.find("Sc") != -1 ):

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.n)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.ellip)
            sherpa.thaw(gcomp.theta)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.ampl)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            sherpa.thaw(gcomp.r0)

            sherpa.fit()
            if (names.display_ds9.find('y') != -1):
                sherpa.image_fit()

            if (psf_gauss.find("Sc") != -1 and psf_gauss.find("P") == -1):
                sherpa.thaw(gcomp.xpos)
                sherpa.thaw(gcomp.ypos)


    sherpa.fit()
    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()

    sherpa.thaw(xrays_mod.ampl)            

    if(names.psf_fudge != -1):
        gpsf.ampl = 1.

    sherpa.fit()
    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()

    sys.stdout.write("\n###########################################\n")
    sherpa.fit(outfile=names.save_path[names.analysis] + "fit_info_"+psf_gauss[:3]+".out", clobber=True)
    print "\n"
    if (names.display_ds9.find('y') != -1):
        sherpa.image_fit()
    print "\n"
    #sherpa.show_model()
    sys.stdout.write("\n###########################################\n")    


    sherpa.show_model()

    # Saving the model 
    elliptical = []
    
    if(psf_gauss.find("G") != -1):
        elliptical   = [gcomp.xpos.val,gcomp.ypos.val,gcomp.fwhm.val,gcomp.ellip.val,gcomp.theta.val,xrays_mod.ampl.val,gcomp.ampl.val]
    if(psf_gauss.find("sh") != -1):    
        elliptical   = [shell.xpos.val,shell.ypos.val,shell.r0.val,shell.width.val,shell.ampl.val,xrays_mod.ampl.val,shell.ampl.val]
    if(psf_gauss.find("dc") != -1):    
        elliptical   = [disc.xpos.val,disc.ypos.val,disc.r0.val,disc.ampl.val,bkg_mod.ampl.val,xrays_mod.ampl.val,disc.ampl.val]
    if(psf_gauss.find("G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1 and psf_gauss.find("Sc") == -1):
        elliptical   = [xrays_mod.ampl.val]
    if(psf_gauss.find("Sc") != -1):
        elliptical   = [gcomp.xpos.val,gcomp.ypos.val,gcomp.r0.val,gcomp.ellip.val,gcomp.theta.val,xrays_mod.ampl.val,gcomp.ampl.val,gcomp.n.val]

    
    sherpa.save_model(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:3] + bkg_option[:3] + psr_alpha + "_model.fits",clobber=True)
    sherpa.save_resid(names.save_path[names.analysis] + names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:3] + bkg_option[:3] + psr_alpha + "_resid.fits",clobber=True)

#==============================================================================
#     sys.stdout.write("\n__________________________________________\nConfidence intervals\nCompute? y/n :\n")
#     confidence   = sys.stdin.readline()
#==============================================================================
    confidence = names.confidence_levels
    if( confidence.find("y") != -1):
        sherpa.conf()
        temporary = sys.stdout
        #sys.stdout = open(names.save_path[names.analysis]+"conf_"+psf_gauss[:3]+".txt", 'w')
        sys.stdout = open(names.save_path[names.analysis] + "conf_" + psf_gauss + ".txt", 'w')
        print(sherpa.get_conf_results())
        sys.stdout.close()
        sys.stdout = temporary       
    sys.stdout.write("\n__________________________________________\n")

    name_file =  names.filename_gamma[names.analysis][:-5] + "_" + psf_gauss[:3] + bkg_option[:3] + psr_alpha

    sys.stdout.write("\n__________________________________________\n")
    
    return psf,psf_gauss,name_file,elliptical

#==============================================================================
# MSH 15-52 with the X-ray templates + factor for given slices           
#==============================================================================
def fit_msh1552_slices(dim,bkg_option,psf_gauss,psr_alpha,boxes,n_slice):

    sherpa.set_method(names.method)
    # Convoluted PSF model
    sherpa.load_psf( "hPSF2", names.path_psf + "rebinned_" + names.filename_psf)
    sherpa.set_psf(hPSF2)
    print "\n - - - - - - - - - - - - \n"
    sherpa.show_psf()
    print "\n - - - - - - - - - - - - \n"
    print(hPSF2)
    #sherpa.image_psf()

    sherpa.set_psf(hPSF2) 
    sherpa.show_kernel()

    #==============================================================================
    #==============================================================================
    
    # X-ray model    
    sherpa.load_table_model( "xrays_mod", names.path_x + "rebinned_400x400.fits" )
    sherpa.show_source()
    print "Xray model :\n",xrays_mod
    
    
    
    # Exposure 
    sherpa.load_table_model("exposure", names.path_expo + names.filename_expo)
    #sherpa.set_exposure(exposure)
    print "Exposure :\n",exposure
    sherpa.show_source()
    #sherpa.set_data("exposure",exposure)
    
    #==============================================================================     
    #==============================================================================

    # Loading the gamma-like candidates image
    sys.stdout.write("\n*** Gamma-like candidates ***\n")
    sherpa.load_image(names.path_template + names.filename_gamma)

    sherpa.set_stat("cash")

    # Using the parameter for the background and adding a gaussian fit for the source

    sherpa.notice2d("box("+str(boxes[0])+","+str(boxes[1])+","+str(boxes[2])+","+str(boxes[3])+","+str(boxes[4])+")")

    psf       = "y"
    psf_gauss = "X" + psf_gauss
    
    if(psf.find("y") != -1 ):
        if(psf_gauss.find("G1") == -1 and psf_gauss.find("2G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1):
            if(bkg_option.find("gam") != -1):
                sherpa.set_full_model(bkg_mod + hPSF2(xrays_mod * exposure))
                
                
    if(bkg_option.find("gam") != -1 ):
        
        bkg_mod.ampl    =   1.1
        sherpa.freeze(bkg_mod.ampl)
        sherpa.freeze(exposure.ampl)
                
        sherpa.fit()
        #sherpa.image_fit()
        sherpa.thaw(bkg_mod.ampl)

        sherpa.fit()
        #sherpa.image_fit()
                                    
        sherpa.thaw(xrays_mod.ampl)   
        #sherpa.thaw(exposure.ampl)
        
        sherpa.fit()
        #sherpa.image_fit()
            
    
    sherpa.thaw(xrays_mod.ampl)            

    sherpa.fit()
    #sherpa.image_fit()
            
    sys.stdout.write("\n###########################################\n")
    sherpa.fit(outfile=names.save_path+"fit_info_slice_"+n_slice+"_"+psf_gauss[:2]+".out", clobber=True)
    print "\n"
    #sherpa.image_fit()
    print "\n"
    sherpa.show_model()
    sys.stdout.write("\n###########################################\n")    
    
    # Saving the model 
    
    elliptical = []
    if(psf_gauss.find("G") == -1 and psf_gauss.find("sh") == -1 and psf_gauss.find("dc") == -1):
        elliptical   = [xrays_mod.ampl.val]

    
    sherpa.save_model(names.save_path + 'sliced/' + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + psr_alpha + "_slice" + n_slice + "_model.fits",clobber=True)
    sherpa.save_resid(names.save_path + 'sliced/' + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + psr_alpha + "_slice" + n_slice + "_resid.fits",clobber=True)


    sys.stdout.write("\n__________________________________________\n")

    name_file =  'sliced/' + names.filename_gamma[:-5] + "_" + psf_gauss[:2] + bkg_option[:2] + psr_alpha + '_slice' + n_slice
    sys.stdout.write("\n__________________________________________\n")
    
    return psf,psf_gauss,name_file,elliptical