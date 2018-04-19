import os, sys
import random, math

import numpy as np
import matplotlib.pyplot as plt


x_plot, y_plot = [], []
x, y = 1., 1.
delta = 0.1
for i in range(100000):
    del_x, del_y   = random.uniform(-delta, delta), random.uniform(-delta, delta)
    if abs(x + del_x) < 1 and abs(y + del_y) < 1:
        x, y = x + del_x, y + del_y
    x_plot.append(x)
    y_plot.append(y)
xyc = range( len( x_plot ) )


plt.clf()
plt.scatter(x_plot,y_plot,c = xyc, marker = '.', s=20, cmap='CMRmap')
plt.axis('equal')

plt.show()
plt.clf()

plt.hist(x_plot,100,normed='True')

plt.show()