file = open("constant.txt", "r")
z = file.readline()

xmin = -20.0
xmax = 20.0
dx = 0.01

xlist = mlab.frange (xmin, xmax, dx)
tlist = mlab.frange (xmin, xmax, dx)

x1 = [2 * np.sqrt(float(z)) * np.cos(t) for t in tlist]
y1 = [np.sqrt(float(z)) * np.sin(t) for t in tlist]

x2 = [1 / 2 * np.sqrt(4 - float(z) ** 2) * np.cos(t) for t in tlist]
y2 = [3 / 2 * np.sqrt(4 - float(z) ** 2) * np.sin(t) for t in tlist]

pylab.plot (x1, y1)
pylab.plot (x2, y2)

plt.show()
