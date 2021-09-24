# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Total ADM mass
M = 1.0
# The mode, as text
mode = "22"
# output data from running merger
data = np.loadtxt("Weyl_integral_" + mode + ".dat")

# make the plot
fig = plt.figure()

# first radius
r0 = 100 #R0 + M*np.log(R0/(2.0*M) - 1.0)
timedata0 = (data[:,0] - r0) / M
fluxdata0 = data[:,1]
plt.plot(timedata0, fluxdata0, '-', lw = 0.75, label="r0")

# second radius
r1 = 150 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata1 = (data[:,0] - r1) / M
fluxdata1 = data[:,3]
plt.plot(timedata1, fluxdata1, '-', lw = 0.75, label="r1")

# third radius
r2 = 200 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata2 = (data[:,0] - r2) / M
fluxdata2 = data[:,5]
plt.plot(timedata2, fluxdata2, '--', lw = 0.75, label="r2")

# make the plot look nice
plt.xlabel("time t / M")
plt.ylabel("Re (Psi4) el, em = " + mode)
plt.xlim(-100, 200)
#plt.ylim(1.85, 1.88)
plt.legend()

# save as png image
filename = "Weyl_" + mode + ".png"
plt.savefig(filename)
