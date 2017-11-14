
# C2EP
#
# Q = C2EP(C) translates the 3x3 direction cosine matrix
# C into the corresponding 4x1 Euler parameter vector Q,
# where the first component of Q is the non-dimensional
# Euler parameter Beta_0 >= 0. Transformation is done
# using the Stanley method.i


import numpy as np


b = np.zeros(1, 4)

def C2EP(input):
    c = input[0:3, 0:3]


    tr = c[1, 1]+c[2, 2]+c[3, 3]
    
    b2 = []
    b2[1] = (1+tr)/4
    b2[2] = (1+2*c[1, 1]-tr)/4
    b2[3] = (1+2*c[2, 2]-tr)/4
    b2[4] = (1+2*c[3, 3]-tr)/4

    [v, i] = np.amax(b2)
    if i == 0:
        def case1():
            b[1] = np.sqrt(b2[1])
            b[2] = (c[2, 3]-c[3, 2])/4/b[1]
            b[3] = (c[3, 1]-c[1, 3])/4/b[1]
            b[4] = (c[1, 2]-c[2, 1])/4/b[1]
            return b

    elif i == 1:
        def case2():
            b[2] = np.sqrt(b2[2])
            b[1] = (c[2, 3]-c[3, 2])/4/b[2]

            if b[1] < 0:
                b[2] = -b[2]
                b[1] = -b[1]

            b[3] = (c[1, 2]+c[2, 1])/4/b[2]
            b[4] = (c[3, 1]+c[1, 3])/4/b[2]
            return b

    elif i == 2:
        def case3():
            b[3] = np.sqrt(b2[3])
            b[1] = (c[3, 1]-c[1, 3])/4/b[3]

            if b[1] < 0:
                b[3] = -b[3]
                b[1] = -b[1]

            b[2] = (c[1, 2]+c[2, 1])/4/b[3]
            b[4] = (c[2, 3]+c[3, 2])/4/b[3]
            return b

    elif i == 3:
        def case4():

            b[4] = np.sqrt(b2[4])
            b[1] = (c[1, 2]-c[2, 1])/4/b[4]

            if b[1] < 0:
                b[4] = -b[4] 
                b[1] = -b[1]

            b[2] = (c[3, 1]+c[1, 3])/4/b[4]
            b[3] = (c[2, 3]+c[3, 2])/4/b[4]
            return b

    bb = b.transpose()
    return bb
