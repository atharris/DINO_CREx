

# C2MRP
# This function translates the 3x3 direction cosine matrix C into the corresponding 3x1 MRP vector S where the MRP
# vector is chosen such that |S| <= 1.
from C2EP import C2EP
import numpy as np

s = np.zeros(1, 3)

def C2MRP(input):

    c = input[0:3, 0:3]

    b = C2EP(c)


    s[1] = b[2]/(1+b[1])
    s[2] = b[3]/(1+b[1])
    s[3] = b[4]/(1+b[1])

    ss = s.transpose()

    return ss


