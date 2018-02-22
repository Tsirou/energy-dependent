import names

import numpy as np
import matplotlib.pyplot as plt

from astropy import wcs
from astropy.io import fits

from conversion import norm_fits

def pulsar_mag(psr_factor):

    # Coordinates of the pulsar
    ra_psr       = names.ra_psr
    dec_psr      = names.dec_psr


    hdu_x        = fits.open(names.path_x + 'rebinned_' + names.dim_for_a + '.fits')
    hdu_gam      = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])


    header      = hdu_gam[names.hdu_ext].header
    w = wcs.WCS(header)

    w.wcs.ctype = ["RA---GLS", "DEC--GLS"]
    coord1 = ra_psr
    coord2 = dec_psr


    x_psr,y_psr  = w.wcs_world2pix(coord1, coord2, names.random_convert_shift)

    # # To check if the coordinate system is adequate
    # print "PSR (pix) : ",x_psr,y_psr
    # plt.clf()
    # plt.imshow(np.log(hdu_x[0].data), cmap='CMRmap')
    #
    # plt.plot(x_psr,y_psr,'+',color='darkgreen',markersize=10,label="PSR B1509-58")
    # plt.show()


    # Normalization  * r^(psr_factor)
    dim_im       = len(hdu_x[0].data)    
    norm_mp      = np.zeros((dim_im,dim_im))
    
    for i in range(0,dim_im):
        for j in range(0,dim_im):      
            norm_mp[i,j]  = hdu_x[0].data[i,j] * (np.sqrt( (j - x_psr)**2 + (i - y_psr)**2))**psr_factor

    hdu_x[0].data              = norm_mp

    negative = names.pulsar_then_rebin
    if(negative.find('n') != -1):
        sum_pixels           = np.sum(hdu_x[0].data, axis=None)
        hdu_x[0].data        = norm_fits(hdu_x[0].data,1./sum_pixels)

    hdu_gam.close()

    hdu_x.writeto(names.path_x + 'rebinned_' + names.dim_for_a + '.fits',overwrite=True)

    hdu_x.close()


    return    