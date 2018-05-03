import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.ndimage as ndimage


from astropy.io import fits
from astropy import wcs
from astropy.modeling.models import Sersic2D
from conversion import factorial,Hermite_polynoms, Sersic_profile, crop
from display import reset_display

reset_display()
fig = plt.figure()
fig.set_size_inches(12.5, 10.5)


# # Hermite polynoms
# degrees  = 5
# x_plus   = (0.01) * np.arange(0,201)
# x_minus  = (-0.01) * np.arange(1,201)
#
# x    = np.concatenate((list(reversed(x_minus)),x_plus),axis=0)
#
#
# for n in range(0,degrees + 1):
#     y   =   Hermite_polynoms(x,n)
#
#     plt.plot(x,y,'-',label="n = "+str(n))
#
#
# plt.legend(fontsize=15, numpoints=1, loc=1)
# plt.show()
#
# plt.clf()


# Sersic parameters
#0.5240
# n         = [0.1, 0.3, 0.5, 0.7, 1.0]#, 2.0, 10.0, 20.0]
# r_0       = 16.74
# epsilon   = 0.3947
# theta     = 3.716
# x_o       = names.X_psr
# y_o       = names.Y_psr
# A         = 0.068
#
#
# # x,y = np.meshgrid(np.arange(400), np.arange(400))
# #
# # mod = Sersic2D(amplitude = A, r_eff = r_0, n=n, x_0=x_o, y_0=y_o, ellip=epsilon, theta=theta)
# # img = mod(x, y)
# # log_img = np.log10(img)
# #
# #
# # plt.xlim(160, 240)
# # plt.ylim(160, 240)
# #
# # plt.imshow(img, origin='lower', cmap="terrain")
# # plt.xlabel('x')
# # plt.ylabel('y')
# # cbar = plt.colorbar()
# # plt.show()
#
#
# for i in n :
#     f   = Sersic_profile(i, r_0, epsilon, theta, x_o, y_o, A)
#
#     plt.plot(f[150:250], label="n = " + str(i)[:4])
#
# plt.legend(fontsize=15, numpoints=1, loc=1)
#
# plt.show()





hdu_gam      = fits.open(names.path_template[names.analysis] + names.filename_gamma[names.analysis])
hdu_exc      = fits.open(names.path_template[names.analysis] + names.filename_excess[names.analysis])

#hdu_bkg_lep  = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis][:-5] + '_smoothed_bkg_'+str(7.0)+'sigma.fits')

hdu_bkg_lep  = hdu_gam[names.hdu_ext].data - hdu_exc[names.hdu_ext].data

hdu_bkg_tmp  = fits.open(names.path_template[names.analysis] + names.filename_bkg[names.analysis])

header = hdu_bkg_tmp[names.hdu_ext].header
w      = wcs.WCS(header)

ax1     = fig.add_subplot(321, projection=w)
ax2     = fig.add_subplot(322, projection=w)
ax3     = fig.add_subplot(323, projection=w)
ax4     = fig.add_subplot(324, projection=w)
ax5     = fig.add_subplot(325, projection=w)
ax6     = fig.add_subplot(326, projection=w)


ax1.imshow(hdu_bkg_lep, cmap="CMRmap")
ax2.imshow(hdu_bkg_tmp[names.hdu_ext].data, cmap="CMRmap")

hdu_lep_tmp   = hdu_bkg_lep / hdu_bkg_tmp[names.hdu_ext].data
hdu_tmp_lep   = hdu_bkg_tmp[names.hdu_ext].data / hdu_bkg_lep


smooth_factor = 7.0
hdu_lep_tmp   = ndimage.gaussian_filter(hdu_bkg_lep, sigma=smooth_factor, order=0) / ndimage.gaussian_filter(hdu_bkg_tmp[names.hdu_ext].data, sigma=smooth_factor, order=0)
hdu_tmp_lep   = ndimage.gaussian_filter(hdu_bkg_tmp[names.hdu_ext].data, sigma=smooth_factor, order=0) / ndimage.gaussian_filter(hdu_bkg_lep, sigma=smooth_factor, order=0)


ax3.imshow(hdu_lep_tmp, cmap="CMRmap")
ax4.imshow(hdu_tmp_lep, cmap="CMRmap")

smooth_factor = 0.0
hdu_lep_tmp   = ndimage.gaussian_filter(hdu_bkg_lep, sigma=smooth_factor, order=0) / ndimage.gaussian_filter(hdu_bkg_tmp[names.hdu_ext].data, sigma=smooth_factor, order=0)
hdu_tmp_lep   = ndimage.gaussian_filter(hdu_bkg_tmp[names.hdu_ext].data, sigma=smooth_factor, order=0) / ndimage.gaussian_filter(hdu_bkg_lep, sigma=smooth_factor, order=0)

ax5.imshow(hdu_lep_tmp, cmap="CMRmap")
ax6.imshow(hdu_tmp_lep, cmap="CMRmap")


plt.show()