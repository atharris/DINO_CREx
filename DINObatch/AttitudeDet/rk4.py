
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
    x0 = x.reshape(6)
    m = m
    sbr = sbr
    wbr = wbr

    m = int(round(m))
    CC = np.zeros(6)

    h = (b - a) / m 
    t = np.zeros([1, m+1])
    t[0, 0] = a
    y = np.column_stack((x0, CC))

    for jj in range(m):
        yjj = y[0:6, jj].transpose()
        k1 = h*feval(yjj, sbr, wbr)

        k2 = h * feval(yjj+k1/2, sbr, wbr)

        k3 = h * feval(yjj+k2/2, sbr, wbr)

        k4 = h * feval(yjj+k3, sbr, wbr)

        y[0:6, jj+1] = yjj + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        s = norm(y[0:3, jj+1])
        if s > 1:
            y[0:3, jj+1] = -y[0:3, jj+1] / (s**2)


    return y[..., jj+1]
