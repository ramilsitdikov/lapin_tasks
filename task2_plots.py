from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy

def drawNorm(Y):
  plt.plot(Y)
  plt.show()

def drawPlot(X, Y, Z, elev=40, azim=-80, minZ = 0.0, maxZ = 0.42):
  fig = plt.figure()
  ax = fig.gca(projection = '3d')
  X, Y = numpy.meshgrid(X, Y)
  surf = ax.plot_surface(X, Y, Z, rstride = 1, cstride = 1, cmap = cm.coolwarm, linewidth = 0, antialiased = False)
  #ax.set_zlim(minZ, maxZ)
  ax.view_init(elev=elev, azim=azim)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  fig.colorbar(surf, shrink = 0.5, aspect = 5)
  plt.show()

h = 0.01

u = [[]]
drawPlot(numpy.arange(0,  1 + h, h), numpy.arange(0,  1 + h, h), u)

p = [[]]
drawPlot(numpy.arange(0,  1 + h, h), numpy.arange(0,  1 + h, h), p, 86, 90)

r = []
drawNorm(r)

iterations = []
drawNorm(iterations)
