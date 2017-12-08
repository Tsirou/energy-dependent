import os
import names

from PyAstronomy import pyasl
from astropy.io import fits



# c          = 2.9992458e8   # m.s-1
# E_nu       = 1.0e12        # eV
#
# wavelength = h * c /E_nu / 1.0e-10
#
# print wavelength
#
# os._exit(0)
#
# hdu_g   = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
#
# image   = hdu_g[names.hdu_ext].data
#
# for i in range(0,len(image)):
#     for j in range(0,len(image)):
#         photons = image[i,j]
#         flux    = pyasl.photons2flux(wavelength,photons)
#
# print(flux)
#
# hdu_g.close()




for comp in range(0,3):
    set_full_model(comp, bkg + psf((gauss2d.g0 + gauss2d.g1 + gauss2d.g2)) * exposure)

    if(comp == 0 ):
        g0.integrate=True
        g1.integrate=False
        g2.integrate=False

    if(comp == 1 ):
        g0.integrate=True
        g1.integrate=True
        g2.integrate=False

    if(comp == 2 ):
        g0.integrate=True
        g1.integrate=True
        g2.integrate=True
