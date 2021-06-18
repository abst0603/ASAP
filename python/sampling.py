import numpy as np
import asap
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

data = np.genfromtxt('synjelly_gudhi.csv', delimiter=',')

samples = asap.sampling(data,radius = 0.15, lr = 1, tau = 3.5)

fig = pyplot.figure()
ax = Axes3D(fig)

ax.scatter(samples[:,0], samples[:,1], samples[:,2])
pyplot.show()