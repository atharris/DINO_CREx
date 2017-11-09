#!/usr/local/bin/python

__author__ = 'Marc Balducci'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################


import sys, os, inspect

import sympy as sym 
import numpy as np
import pdb
################################################################################
#                  S E C O N D A R Y     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

################################################################################
#             U N I T     T E S T     C A S E     F U N C T I O N:
################################################################################

# -------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------

################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    # define all symbols necessary
    # x position of spacecraft
    x_sc = sym.Symbol('x_sc')
    # y position of spacecraft
    y_sc = sym.Symbol('y_sc')
    # z position of spacecraft
    z_sc = sym.Symbol('z_sc')
    # x position of object
    x_ob = sym.Symbol('x_ob')
    # y position of object
    y_ob = sym.Symbol('y_ob')
    # z position of object
    z_ob = sym.Symbol('z_ob')
    # focal length (mm)
    FoL  = sym.Symbol('FoL')
    # right ascension 
    RA   = sym.Symbol('RA')
    # declination
    dec  = sym.Symbol('dec')
    # twist
    twist = sym.Symbol('twist')
    # the x dimension conversion of a pixel
    Kx   = sym.Symbol('Kx')
    # the y dimension conversion of a pixel
    Ky   = sym.Symbol('Ky')
    # x to pixel direction
    Dx   = sym.Symbol('Dx')
    # y to pixel direction
    Dy   = sym.Symbol('Dy')
    # pixel center offset
    p0   = sym.Symbol('p0')
    # line center offset
    l0   = sym.Symbol('l0')

    # collect syms into position vectors 
    state_spacecraft     = sym.Matrix([[x_sc],[y_sc],[z_sc]])
    state_object         = sym.Matrix([[x_ob],[y_ob],[z_ob]])

    # find distance vector (from spacecraft to object)
    spacecraft_to_object = state_object - state_spacecraft

    # rotation matrix from inertial (I) to camera (TV)
    T_TV_I = rot3( RA ) * rot2( dec ) * rot1( twist )

    # unit vector for pointing in inertial
    A_hat_I = spacecraft_to_object / sym.sqrt( spacecraft_to_object[0]**2 + \
              spacecraft_to_object[1]**2 + spacecraft_to_object[2]**2 )

    # rotate unit vector from inertial to camera
    A_hat_TV = T_TV_I * A_hat_I

    # pixel conversion
    p = Kx * Dx * FoL / A_hat_TV[2] * A_hat_TV[0] + p0
    # line conversion
    l = Ky * Dy * FoL / A_hat_TV[2] * A_hat_TV[1] + l0

    # differentiate pixel and line wrt to spacecraft position 
    p_deriv_x = sym.diff( p, x_sc )
    p_deriv_y = sym.diff( p, y_sc )
    p_deriv_z = sym.diff( p, z_sc )
    
    l_deriv_x = sym.diff( l, x_sc )
    l_deriv_y = sym.diff( l, y_sc )
    l_deriv_z = sym.diff( l, z_sc )

    # turn symbolics into functions
    p_deriv_x_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), p_deriv_x )

    p_deriv_y_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), p_deriv_y )

    p_deriv_z_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), p_deriv_z )

    l_deriv_x_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), l_deriv_x )

    l_deriv_y_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), l_deriv_y )

    l_deriv_z_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), l_deriv_z )

    A_hat_I_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), A_hat_I )

    A_hat_TV_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), A_hat_TV )

    p_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), p )

    l_function = sym.lambdify( ( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 ), l )

    # redefine variables as floats
    x_sc  = 1000.
    y_sc  = 0.
    z_sc  = 0.
    x_ob  = 1200.
    y_ob  = 1000.
    z_ob  = 450.
    FoL   = 100.
    RA    =   0.
    dec   =   0.
    twist = np.pi/2.
    Kx    = 1./2.5
    Ky    = 1./10.
    Dx    = 1.
    Dy    = 1.
    p0    = 2048./2.
    l0    = 512./2.

    # evaluate all at specified values
    p_deriv_x_test = p_deriv_x_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    p_deriv_y_test = p_deriv_y_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    p_deriv_z_test = p_deriv_z_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    l_deriv_x_test = l_deriv_x_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    l_deriv_y_test = l_deriv_y_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    l_deriv_z_test = l_deriv_z_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    A_hat_I_test = A_hat_I_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    A_hat_TV_test = A_hat_TV_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0 )

    p_test = p_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0  )
    
    l_test = l_function( x_sc, y_sc, z_sc, x_ob, y_ob,\
      z_ob, FoL, RA, dec, twist, Kx, Ky, Dx, Dy, p0, l0  )

    return [ p_test, l_test, p_deriv_x_test, p_deriv_y_test, p_deriv_z_test, \
             l_deriv_x_test, l_deriv_y_test, l_deriv_z_test ]

def rot1( angle ):
    rot1_mat = sym.Matrix( [[1,0,0],[0,sym.cos(angle),-sym.sin(angle)],[0,sym.sin(angle),sym.cos(angle)]] )
    return rot1_mat
def rot2( angle ):
    rot2_mat = sym.Matrix( [[sym.cos(angle),0,-sym.sin(angle)],[0,1,0],[sym.sin(angle),0,sym.cos(angle)]] )
    return rot2_mat
def rot3( angle ):
    rot3_mat = sym.Matrix( [[sym.cos(angle),sym.sin(angle),0],[-sym.sin(angle),sym.cos(angle),0],[0,0,1]] )
    return rot3_mat

if __name__ == "__main__":
    main()
