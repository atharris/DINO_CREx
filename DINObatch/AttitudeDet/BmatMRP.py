
# BmatMRP(Q)
#
#	B = BmatMRP(Q) returns the 3x3 matrix which relates the 
#	body angular velocity vector w to the derivative of
#	MRP vector Q.  
#	
#		dQ/dt = 1/4 [B(Q)] w
#	



import numpy as np

b = np.zeros(3, 3)

def BmatMRP(input):
    
    s = input[0:3]

    # B = BmatMRP(s) returns the 3x3 matrix which relates the 
    # body angular velocity vector w to the derivative of
    # MRP vector s.  
    # ds/dt = 1/4 [b[s)] w
    #
    

    s2 = s.transpose()*s
    b[1, 1] = 1-s2+2*s[1]*s[1] 
    b[1, 2] = 2*(s[1]*s[2]-s[3]) 
    b[1, 3] = 2*(s[1]*s[3]+s[2]) 
    b[2, 1] = 2*(s[2]*s[1]+s[3]) 
    b[2, 2] = 1-s2+2*s[2]*s[2] 
    b[2, 3] = 2*(s[2]*s[3]-s[1]) 
    b[3, 1] = 2*(s[3]*s[1]-s[2]) 
    b[3, 2] = 2*(s[3]*s[2]+s[1]) 
    b[3, 3] = 1-s2+2*s[3]*s[3] 
    
    return b 