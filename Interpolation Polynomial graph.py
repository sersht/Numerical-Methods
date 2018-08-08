# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import scipy.interpolate as sp
from numpy import array
import pylab
from matplotlib import mlab
from scipy import interpolate

# coordinates on x-axis
with open("points_for_interpolation.pr4") as fx:
    lines = fx.read().splitlines()
points = [float(item) for item in lines]

# Chebyshev points and values of polynomial there
with open("cheb_y.pr4") as fy:
    lines = fy.read().splitlines()
chebY = [float(item) for item in lines]

# Equidistant points and values of polynomial there
# равноудаленные точки и значения полинома в них
with open("equidist_y.pr4") as fy:
    lines = fy.read().splitlines()
equidistY = [float(item) for item in lines]

# Function graph
xmin = -100.0
xmax = 100.0
dx = 0.001

xlist = mlab.frange (xmin, xmax, dx)
ylist = [np.exp(-z * z) for z in xlist]
plt.plot(xlist, ylist, 'r')

text = 'quadratic'

# Chebyshev graph
chebXarray = array(points)
chebYarray = array(chebY)
f = interpolate.interp1d(chebXarray, chebYarray, kind=text)
xnew = np.linspace(chebXarray.min(), chebXarray.max())

# Equidistant graph
x_ = array(points)
y_ = array(equidistY)
F = interpolate.interp1d(x_, y_, kind=text)
Xnew = np.linspace(x_.min(), x_.max())

plt.plot(xnew, f(xnew), 'g')
plt.plot(Xnew, F(Xnew), 'b')

pylab.legend ( ("f(x)", "Chebyshev", "Equidistant") )

plt.axis([-2, 2, 0, 2])
plt.show()
