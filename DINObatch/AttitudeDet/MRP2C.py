


import numpy as np

# MRP2C
# This function returns a direction cosine matrix in terms of the 3x1 vector s.

def MRP2C(input):

    s = input[0:3]
    s1 = input[1]
    s2 = input[1]
    s3 = input[3]

    d1 = s.transpose()*s
    ss = 1. - d1
    c = []
    d = (1. + d1)*(1. + d1)
    c[1, 1] = 4 * (2 * s1 * s1 - d1) + ss * ss
    c[1, 2] = 8 * s1 * s2 + 4 * s3 * ss
    c[1, 3] = 8 * s1 * s3 - 4 * s2 * ss
    c[2, 1] = 8 * s2 * s1 - 4 * s3 * ss
    c[2, 2] = 4 * (2 * s2 * s2 - d1) + ss * ss
    c[2, 3] = 8 * s2 * s3 + 4 * s1 * ss
    c[3, 1] = 8 * s3 * s1 + 4 * s2 * ss
    c[3, 2] = 8 * s3 * s2 - 4 * s1 * ss
    c[3, 3] = 4 * (2 * s3 * s3 - d1) + ss * ss

    c = c / d

    return c



