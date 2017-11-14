
# #---------------------------------------------------------------------------
# RK4   Runge-Kutta solution for y' = f(t,y) with y(0) = X0. (a = 0)
# Sample call
#   [T,Y] = rk4('f',a,b,X0,m,options)
# Inputs
#   feval    name of the function
#   a    left  endpoint of [a,b]
#   b    right endpoint of [a,b]
#   X0   initial value
#   m    number of steps
#   options    other things passed through in a struct
# Return
#   T    solution: vector of abscissas
#   Y    solution: vector of ordinates
#--------------------------------------------------------------------------


import numpy as np
from scfun import scfun
import pdb

def norm( input ):
  norm = np.sqrt(sum(np.square(input)))
  return norm


feval = scfun


def rk4(a, bb, x, m, sbr, wbr):
    a = a
    b = bb
    x0 = x.reshape(6, 1)
    m = m
    sbr = sbr
    wbr = wbr

    m = int(round(m))

    h = (b - a) / m 
    t = np.zeros([1, m+1])
    t[0, 0] = a
    y = np.array([x0[0:6, 1], [np.zeros([6, 1])]])
    pdb.set_trace()
    for j in range(m):

        yj = y[0:5, j-1].transpose()

        k1 = h*feval(yj, sbr, wbr)

        k2 = h * feval(yj+k1/2, sbr, wbr)

        k3 = h * feval(yj+k2/2, sbr, wbr)

        k4 = h * feval(yj+k3, sbr, wbr)
        y[j] = yj + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        s = norm(y[0:3, j])
        if s > 1:
            y[0:3, j] = -y[0:3, j] / (s**2)

    return y
