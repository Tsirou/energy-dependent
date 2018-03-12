import names

import numpy as np
import matplotlib.pyplot as plt
from conversion import rebin_med
from conversion import fill_fits
from conversion import norm_fits

from astropy import wcs
from astropy.io import fits


def xray_map_resize(negative):

    if(names.pulsar_then_rebin.find('y') != -1):
        hdu_x                      = fits.open(names.path_x + 'rebinned_' + names.dim_for_a + '.fits')
    else:
        hdu_x                      = fits.open(names.path_x + names.filename_x[names.analysis])

    im                         = hdu_x[0].data

    # Rebinning by a factor of 9 --> 100x100
    im                         = rebin_med( im, 9 )    
    # Filling a 400x400 grid
    im                         = fill_fits(im,400)


    negative = names.pulsar_then_rebin
    if(negative.find('y') != -1):
        sum_pixels           = np.sum(im, axis=None)
        im                   = norm_fits(im,1./sum_pixels)

    hdu_x[0].data              = im


    hdu_x.writeto(names.path_x + 'rebinned_'+ names.dim_for_a +'.fits',clobber=True)

    hdu_x.close()

    return


def xray_header():

    hdu_xray                   = fits.open(names.path_x + 'rebinned_' + names.dim_for_a + '.fits')
    hdu_x                      = fits.open(names.path_x + names.filename_x[names.analysis])
    hdu                        = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])

    dim_init                   = len(hdu_x[0].data)

    an = 1

    # FITS HEADER switches
    hdu_xray[0].header['crpix1']  = len(hdu[names.hdu_ext].data)/2. + 0.5 + (hdu_x[0].header['crpix1'] - (dim_init/2. + 0.5))/(9 * an)
    hdu_xray[0].header['crpix2']  = len(hdu[names.hdu_ext].data)/2. + 0.5 + (hdu_x[0].header['crpix2'] - (dim_init/2. + 0.5))/(9 * an)
    hdu_xray[0].header['cdelt1']  = hdu[names.hdu_ext].header['cdelt1']
    hdu_xray[0].header['cdelt2']  = hdu[names.hdu_ext].header['cdelt2']
#===========S===================================================================

    hdu_x.close()
    hdu.close()

    hdu_xray.writeto(names.path_x + 'rebinned_' + names.dim_for_a + '.fits',overwrite=True)

    hdu_xray.close()

    return

