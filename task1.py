import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator,FormatStrFormatter
import matplotlib.pyplot as plt

def zadachaOPrepyatstvii(h,f,fi,xi,sigma,epsilon=0.0000001):
  m = int(1/h)
  h2 = h * h
  u = numpy.zeros(shape=(m+1, m+1))
  gamma = numpy.zeros(shape=(m+1, m+1))
  while True:
    for i in xrange(1, m):
      for j in xrange(1, m):
        x = (i * h, j * h)
        v = (1 - sigma) * u[i][j] + sigma * 0.25 * (u[i - 1][j]+ u[i + 1][j] + u[i][j - 1] + u[i][j + 1] + h2 * f(i * h,j * h))
        fi_x = fi(x)
        if v < fi_x:
          u[i][j] =fi_x
        else:
          xi_x =xi(x)
        u[i][j] =xi_x if v >xi_x else v
        gamma[i][j] = 4/ (h2 * sigma)* (v - u[i][j])
        norm_r = 0
    for i in xrange(1,m):
      for j in xrange(1,m):
        r = 1/ h2 * (4* u[i][j]- u[i-1][j]- u[i+1][j]- u[i][j-1]- u[i][j+1])+ gamma[i][j]- f(i*h,j*h)
        norm_r +=h2 * r * r
        print pow(norm_r, 0.5)
    if pow(norm_r, 0.5) <epsilon:
      break
  return u

def drawPlot(X,Y,Z,minZ=-1.01,maxZ=2.01):
  fig =plt.figure()
  ax =fig.gca(projection='3d')
  X,Y =numpy.meshgrid(X,Y)
  surf =ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
  ax.set_zlim(minZ,maxZ)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  fig.colorbar(surf,shrink=0.5,aspect=5)
  plt.show()

def f(x,y):
  if x > 0.3 and x < 0.7 and y > 0.3 and y < 0.7:
    return -60
  return 40

h = 0.02
res =zadachaOPrepyatstvii(h,f, lambda x: 0, lambda x: 10000000.7, 1.9, 0.02)
drawPlot(numpy.arange(0, 1+ h,h),numpy.arange(0, 1+ h,h),res)
