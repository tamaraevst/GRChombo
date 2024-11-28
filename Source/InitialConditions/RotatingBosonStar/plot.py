import numpy as np
from matplotlib import pyplot, rc

data = np.loadtxt("new-dthomega.dat")

shw = pyplot.imshow(data, origin="lower", extent=(0,1,0,1), aspect="auto")
bar = pyplot.colorbar(shw)
pyplot.xlabel(r"$x=r/(1-r)$")
pyplot.ylabel(r"$\theta/\pi$")
pyplot.show()
