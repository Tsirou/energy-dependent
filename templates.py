import os
import names

import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt

from conversion import fill_fits,crop
from conversion import rebin_red,rebin_med,rebin_med_non_avg

from astropy import wcs
from astropy.io import fits


# Template background by subtracting the excess map to the gamma-like candidates one
def leptonic_bkg(smooth_factor):
    hdu_gam_temp   = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
    hdu_exc_temp   = fits.open(names.path_template[names.analysis] + names.filename_excess[names.analysis])

    hdu_background = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis])


    header         = hdu_background[1].header
    w = wcs.WCS(header)
    w.wcs.ctype = ["RA---GLS", "DEC--GLS"]

    excess         = fill_fits(hdu_exc_temp[1].data, len(hdu_gam_temp[1].data))
    #excess         = hdu_exc_temp[1].data
    bkgg           = hdu_gam_temp[1].data - excess

    hdu_background[1].data = bkgg
    hdu_background.writeto(names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_unsmoothed_bkg.fits',clobber=True)


    img      = ndimage.gaussian_filter(bkgg, sigma=smooth_factor, order=0)
    hdu_background[1].data = img
    hdu_background.writeto(names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_smoothed_bkg_'+str(smooth_factor)+'sigma.fits',clobber=True)


    hdu_gam_temp.close()
    hdu_exc_temp.close()
    hdu_background.close()                       


    return


# Cannot remember what this was for...
def resolution_maps(fitsname,factor_fits):

    for name in fitsname:
        # should be : bkg ; gamma ; exposure

        if(name.find('exposure') != -1):
            temporary     = names.hdu_ext
            names.hdu_ext = temporary - 1

        hdu_t      = fits.open(name)

        if(name.find('exposure') != -1):
            hdu_t[names.hdu_ext].data = rebin_med(hdu_t[names.hdu_ext].data,factor_fits)
        else:
            hdu_t[names.hdu_ext].data = rebin_med_non_avg(hdu_t[names.hdu_ext].data, factor_fits)

        # FITS HEADER switches
        hdu_t[names.hdu_ext].header['crpix1'] = (len(hdu_t[names.hdu_ext].data) / 2. / factor_fits) + 0.5 + (hdu_t[names.hdu_ext].header['crpix1'] - (len(hdu_t[names.hdu_ext].data) / 2. + 0.5)) / (factor_fits)
        hdu_t[names.hdu_ext].header['crpix2'] = (len(hdu_t[names.hdu_ext].data) / 2. / factor_fits) + 0.5 + (hdu_t[names.hdu_ext].header['crpix2'] - (len(hdu_t[names.hdu_ext].data) / 2. + 0.5)) / (factor_fits)

        hdu_t[names.hdu_ext].header['cdelt1'] = hdu_t[names.hdu_ext].header['cdelt1'] * factor_fits
        hdu_t[names.hdu_ext].header['cdelt2'] = hdu_t[names.hdu_ext].header['cdelt2'] * factor_fits

        hdu_t.writeto(name[:-5] + '_rebinned.fits', overwrite=True)
        hdu_t.close()

        if(name.find('exposure') != -1):
            names.hdu_ext = temporary


    return


# Counting the counts within a given radius
def photon_count(radius):

    photons = 0

    hdu_p     = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
    #hdu_b     = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_smoothed_bkg_'+str(7.0)+'sigma.fits')
    hdu_b = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_unsmoothed_bkg.fits')
    #hdu_b = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis])

    hdu_e = fits.open(names.path_template[names.analysis] + names.filename_excess[names.analysis])

    binsize   = hdu_p[names.hdu_ext].header['cdelt2']

    radius    = int(radius / binsize) # radius in pixel
    print "\nRadius of selection (in pixels) : ",radius

    subtract  = hdu_p[names.hdu_ext].data
    selection = crop(subtract, radius*2)

    plt.imshow(selection,cmap='CMRmap')
    plt.show()

    circle    = np.zeros((len(selection),len(selection)))
    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if( (i)**2 + (j)**2 <= radius**2):
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 - j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 - j] = 1

    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if(circle[i,j] == 0):
                selection[i,j] = 0

    plt.imshow(selection,cmap='CMRmap')
    plt.show()

    photons = np.sum(selection, axis=None)

    print "\nPhotons based on only the raw count map :",photons

    photons   = 0
    subtract  = hdu_p[names.hdu_ext].data[:,:] - hdu_b[names.hdu_ext].data[:,:]
    selection = crop(subtract, radius*2)

    circle    = np.zeros((len(selection),len(selection)))
    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if( (i)**2 + (j)**2 <= radius**2):
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 - j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 - j] = 1

    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if(circle[i,j] == 0):
                selection[i,j] = 0


    photons   = np.sum(selection,axis=None)
    print "\nPhotons based on raw count map - bkg :",photons

    photons   = 0
    selection = crop(hdu_e[names.hdu_ext].data,radius*2)

    circle    = np.zeros((len(selection),len(selection)))
    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if( (i)**2 + (j)**2 <= radius**2):
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 + j] = 1
                circle[len(selection)/2 - 1 + i, len(selection)/2 - 1 - j] = 1
                circle[len(selection)/2 - 1 - i, len(selection)/2 - 1 - j] = 1

    for i in range(0,len(selection)):
        for j in range(0,len(selection)):
            if(circle[i,j] == 0):
                selection[i,j] = 0


    plt.imshow(selection,cmap='CMRmap')
    plt.show()


    photons   = np.sum(selection,axis=None)

    print "\nPhotons based on the excess template :",photons

    hdu_p.close()
    hdu_b.close()
    hdu_e.close()

    return photons