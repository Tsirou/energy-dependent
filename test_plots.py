import names

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


from astropy.modeling.models import Sersic2D
from conversion import factorial,Hermite_polynoms, Sersic_profile
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

n         = 0.5240
r_0       = 16.74
epsilon   = 0.3947
theta     = 3.716
x_o       = names.X_psr
y_o       = names.Y_psr
A         = 0.068


# x,y = np.meshgrid(np.arange(400), np.arange(400))
#
# mod = Sersic2D(amplitude = A, r_eff = r_0, n=n, x_0=x_o, y_0=y_o, ellip=epsilon, theta=theta)
# img = mod(x, y)
# log_img = np.log10(img)
#
#
# plt.xlim(160, 240)
# plt.ylim(160, 240)
#
# plt.imshow(img, origin='lower', cmap="terrain")
# plt.xlabel('x')
# plt.ylabel('y')
# cbar = plt.colorbar()
# plt.show()



f   = Sersic_profile(n, r_0, epsilon, theta, x_o, y_o, A)

plt.plot(f)
plt.show()

