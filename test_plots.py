import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from conversion import factorial,Hermite_polynoms


degrees  = 5
x_plus   = (0.01) * np.arange(0,201)
x_minus  = (-0.01) * np.arange(1,201)

x    = np.concatenate((list(reversed(x_minus)),x_plus),axis=0)


for n in range(0,degrees + 1):
    y   =   Hermite_polynoms(x,n)

    plt.plot(x,y,'-',label="n = "+str(n))


plt.legend(fontsize=15, numpoints=1, loc=1)
plt.show()