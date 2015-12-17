import numpy
from decimal import *
import math

def create_zero_matrix(dimension):
  return numpy.zeros(shape = (dimension, dimension))

def show_matrix(matrix):
  for row in xrange(len(matrix)):
    for col in xrange(len(matrix[row])):
      print '[{0}, {1}] = {2}'.format(row, col, matrix[row][col])

def some_noname_function(row, col, step):
  x = Decimal(row * step - 0.5)
  y = Decimal(col * step - 0.5)
  if x * x + y * y <= 0.1:
    return 0
  return 4

def solve_problem(step, sigma, epsilon, tau):
  M = int(round(Decimal(1/step)))
  N = M + 1

  step_square = Decimal(M * M)
  u = numpy.zeros(shape = (N, N))
  F = numpy.zeros(shape = (N, N))
  p1 = numpy.zeros(shape = (N, N))
  p2 = numpy.zeros(shape = (N, N))
  l1 = numpy.zeros(shape = (N, N))
  l2 = numpy.zeros(shape = (N, N))

  #outer while loop
  out_iterator = 0
  while  True:

    for row in xrange(0, N):
      for col in xrange(0, N):
        module = math.sqrt(pow(l1[row][col], 2) + pow(l2[row][col], 2))
        print module
        if module <= 1:
          p1[row][col] = (l1[row][col]/module)
          p2[row][col] = (l2[row][col]/module)
        else:
          p1[row][col] = (l1[row][col]/module)
          p2[row][col] = (l2[row][col]/module)

    # print "Matrix p1:"
    # show_matrix(p1)
    # print "Matrix p2:"
    # show_matrix(p2)

    for row in xrange(1, M):
      for col in xrange(1, M):
        Su = (2 * (u[row][col] -u[row - 1][col - 1]))/step
        Lp = 2 * (p1[row][col] - p1[row + 1][col] + p2[row][col] - p2[row][col + 1])/step - (l1[row][col] - l1[row + 1] + l2[row][col] - l2[row][col + 1])/step
        temp = Su/2 + Lp/2 + some_noname_function(row, col, step)/2
        print temp
        F[row][col] = temp

    #inner while loop
    in_iterator = 0
    while True:

      for row in xrange(1, M):
        for col in xrange(1, M):
          u[row][col] = (1 - sigma) * u[row][col] + sigma / 4 * (u[row - 1][col] + u[row + 1][col]  + u[row][col - 1]  + u[row ][col + 1] + step_square * F[row][col])

      norm_r = 0
      for row in xrange (1, M):
        for col in xrange (1, M):
          r = Decimal(4 * u[row][col] - u[row - 1][col] - u[row + 1][col]  - u[row][col - 1]  - u[row ][col + 1]) / step_square - F[row][col]
          norm_r += step_square * pow(r, 2)

      if math.sqrt(norm_r) < epsilon:
        print ", inner iterator = ", in_iterator
        break
    #end of inner while loop

    for row in xrange(1, N):
      for col in xrange (1, N):
        l1[row][col] = l1[row][col] - tau * (p1[row][col] - (u[row][col] - u[row - 1][col]) / step)
        l2[row][col] = l2[row][col] - tau * (p2[row][col] - (u[row][col] - u[row][col - 1]) / step)

    for row in xrange(1, N):
      for col in xrange (1, N):
        d1 = Decimal(p1[row][col] - (u[row][col] - u[row - 1][col] / step))
        d2 = Decimal(p2[row][col] - (u[row][col] - u[row][col - 1] / step))
        r += pow(d1, 2) + pow(d2, 2)

    r = step * math.sqrt(r)
    print ", r = ", r
    if r < epsilon:
      break
  #end of outer while loop

  print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in u]))

  for row in xrange(1, N):
    for col in xrange(1, N):
      print math.sqrt(pow(p1[row][col], 2) + pow(p2[row][col], 2))

  return out_iterator


iterations_number = solve_problem(0.01, 1.92, 0.01, 1.9)
