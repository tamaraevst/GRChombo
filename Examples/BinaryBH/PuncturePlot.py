import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

L = 256
r = 6

fig = plt.figure()
data = np.loadtxt("BBH2_chkExtraction2.txt")
x1 = data[:,1]-L
y1 = data[:,2]-L
plt.plot(x1,y1)
x2 = data[:,4]-L
y2 = data[:,5]-L
plt.plot(x2,y2)
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.ylim(-r-1,r+1)
plt.savefig("PunctureTracks.pdf")
