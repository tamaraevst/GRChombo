# A simple python script to plot the puncture
# tracks over time. This helps with setting up
# circular orbits. The params.txt file should
# give around 8-9 orbits before merger.

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# location of center of grid
L = 256
# half the separation of punctures
r = 6
# output data from running merger
data = np.loadtxt("../BinaryBHChk_Punctures.dat")

# make the plot
fig = plt.figure()

# first puncture
x1 = data[:,1]-L
y1 = data[:,2]-L
plt.plot(x1,y1)

# second puncture
x2 = data[:,4]-L
y2 = data[:,5]-L
plt.plot(x2,y2)

# make the plot look nice
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.ylim(-r-1,r+1)

# save as png image
plt.savefig("PunctureTracks.png")
