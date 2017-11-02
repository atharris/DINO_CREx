import numpy as np
import math

##################################################
##################################################


# theta in radians
def eulerDCM_321(theta1,theta2,theta3):
    'Input theta1, theta2, theta3 (rad) to get DCM for Euler 3-2-1 Rotation'
    c1 = math.cos(theta1)
    s1 = math.sin(theta1)
    c2 = math.cos(theta2)
    s2 = math.sin(theta2)
    c3 = math.cos(theta3)
    s3 = math.sin(theta3)
    C = np.array([[c2*c1, c2*s1, -s2],
                  [s3*s2*c1-c3*s1, s3*s2*s1+c3*c1, s3*c2],
                  [c3*s2*c1+s3*s1, c3*s2*s1-s3*c1, c3*c2]])
    return C


##################################################
##################################################

#theta in radians
def eulerDCM_313(theta1,theta2,theta3):
    'Input theta1, theta2, theta3 (rad) to get DCM for Euler 3-1-3 Rotation'
    c1 = math.cos(theta1)
    s1 = math.sin(theta1)
    c2 = math.cos(theta2)
    s2 = math.sin(theta2)
    c3 = math.cos(theta3)
    s3 = math.sin(theta3)
    C = np.array([[c3*c1-s3*c2*s1, c3*s1+s3*c2*c1, s3*s2],
                 [-s3*c1-c3*c2*s1, -s3*s1+c3*c2*c1, c3*s2],
                 [s2*s1, -s2*c1, c2]])
    return C


##################################################
##################################################
# Converts right ascension and declination to unit vector

# Inputs:   radec   tuple of right ascension and declination [deg]
# Outputs   ehat    1x3 array of unit vector

def radec_to_unitvector(radec):

    ehat = np.array([   math.cos(math.radians(radec[1])) * math.cos(math.radians(radec[0])),
                        math.cos(math.radians(radec[1])) * math.sin(math.radians(radec[0])),
                        math.sin(math.radians(radec[1]))
                        ])

    return ehat


##################################################
##################################################
# Converts right ascension and declination to unit vector

# Inputs    ehat    1x3 array of unit vector
# Outputs:  radec   tuple of right ascension and declination [deg]
#                   right ascensions is 0 < RA < 360

def unitvector_to_radec(ehat):

    e1 = ehat[0]
    e2 = ehat[1]
    e3 = ehat[2]

    hplane = math.sqrt(e1**2 + e2**2)

    dec = math.degrees(math.atan2(e3, hplane))
    ra = math.degrees(math.atan2(e2, e1))

    if ra < 0:
        ra = ra + 360

    radec = (ra, dec)

    return radec