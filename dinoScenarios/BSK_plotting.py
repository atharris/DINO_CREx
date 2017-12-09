import sys, os, inspect
import matplotlib.pyplot as plt
from os.path import expanduser
import numpy as np
from numpy import linalg as la

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')


try:
    import macros as mc
except ImportError:
    import Basilisk.utilities.macros as mc

import socket

home = expanduser("~")
try:
    hostIP = socket.gethostbyname(socket.gethostname())
    print 'Host Name: ', home
    print 'Host IP: ',hostIP
except:
    hostIP = None


m2km = 1.0/1000.0

# ------------------------------------------------- PLOT SAVING ------------------------------------------------- #

paperPath = home + '/Desktop/ISSFD_paper/Figures/SingleProcess' # define path to figures folder
arePlotsSaved = False
def save_plot(title):
    #plt.rcParams['figure.figsize'] = 3.5, 3.0
    plt.rcParams.update({'font.size': 7})
    plt.savefig(paperPath + "/" + title + ".pdf", bbox_inches='tight')

# --------------------------------- COMPONENTS & SUBPLOT HANDLING ----------------------------------------------- #
color_x = 'dodgerblue'
color_y = 'salmon'
color_z = 'lightgreen'

def plot3components(vec):
    time = vec[:, 0] * mc.NANO2MIN
    plt.xlabel('Time, min')
    plt.plot(time, vec[:, 1], color_x)
    plt.plot(time, vec[:, 2], color_y)
    plt.plot(time, vec[:, 3], color_z)

def plot3Covar(covar):
    time = covar[:, 0] * mc.NANO2MIN
    plt.xlabel('Time, min')
    plt.plot(time, 3*np.sqrt(covar[:, 1]), 'b--')
    plt.plot(time, 3*np.sqrt(covar[:, 2]), 'r--')
    plt.plot(time, 3*np.sqrt(covar[:, 3]), 'g--')
    plt.plot(time, -3*np.sqrt(covar[:, 1]), 'b--')
    plt.plot(time, -3*np.sqrt(covar[:, 2]), 'r--')
    plt.plot(time, -3*np.sqrt(covar[:, 3]), 'g--')

def plot_postfit(postFits, noise):
    time = postFits[:, 0] * mc.NANO2MIN
    plt.plot(time, postFits[:,1], 'b.')
    plt.plot(time, 3*noise*np.ones(len(postFits[:, 0])), 'r--')
    plt.plot(time, -3*noise*np.ones(len(postFits[:, 0])), 'r--')
    plt.xlabel('Time, min')
    plt.ylabel('PostFit')

def plot_sigma(sigma):
    plot3components(sigma)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.ylabel('MRP')

def plot_omega(omega):
    plot3components(omega)
    plt.ylabel('Angular Rate, rad/s')
    plt.legend(['$\omega_1$', '$\omega_2$', '$\omega_3$'])

def plot_sigmaCovar(sigma_err, covar):
    plot3components(sigma_err)
    plot3Covar(covar)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.ylabel('MRP')

def plot_omegaCovar(omega_err, covar):
    plot3components(omega_err)
    plot3Covar(covar)
    plt.ylabel('Angular Rate, rad/s')
    plt.legend(['$\omega_1$', '$\omega_2$', '$\omega_3$'])

def subplot_sigma(subplot, sigma):
    plot3components(sigma)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.ylabel('MRP')

def subplot_omega(subplot, omega):
    plot3components(omega)
    plt.ylabel('Angular Rate, rad/s')
    plt.legend(['$\omega_1$', '$\omega_2$', '$\omega_3$'])


# ------------------------------------- MAIN PLOT HANDLING ------------------------------------------------------ #

def plot_bodyTorque(Lb):
    plt.figure()
    plot3components(Lb)
    plt.title('Body Torque $L_b$')
    plt.ylabel('Body Torque, $N \cdot m$')
    plt.legend(['$L_{b,1}$', '$L_{b,2}$', '$L_{b,3}$'])

def plot_controlTorque(Lr):
    plt.figure()
    plot3components(Lr)
    plt.ylabel('Control Torque, $N \cdot m$')
    plt.legend(['$L_{r,1}$', '$L_{r,2}$', '$L_{r,3}$'])
    if arePlotsSaved:
        save_plot("Lr")
    else:
        plt.title('Control Torque $L_r$')
    return

def plot_trackingError(sigma_BR, omega_BR_B):
    if arePlotsSaved:
        plt.figure()
        plot_sigma(sigma_BR)
        save_plot("sigma_BR")
        plt.figure()
        plot_omega(omega_BR_B)
        save_plot("omega_BR_B")
    else:
        plt.figure()
        plt.subplot(211)
        plt.title('Att Error: $\sigma_{BR}$')
        plot_sigma(sigma_BR)
        plt.subplot(212)
        plt.title('Rate Error: $^B{\omega_{BR}}$')
        plot_omega(omega_BR_B)
    return

def plot_attitudeGuidance(sigma_RN, omega_RN_N):
    if arePlotsSaved:
        plt.figure()
        plot_sigma(sigma_RN)
        plt.ylim([-1.0, 1.0])
        save_plot("sigma_RN")
        plt.figure()
        plot_omega(omega_RN_N)
        save_plot("omega_RN_N")
    else:
        plt.figure()
        plt.subplot(211)
        plt.title('Ref Att: $\sigma_{RN}$')
        plot_sigma(sigma_RN)
        plt.ylim([-1.0, 1.0])
        plt.subplot(212)
        plt.title('Ref Rate: $^N{\omega_{RN}}$')
        plot_omega(omega_RN_N)
    return

def plot_rotationalNav(sigma_BN, omega_BN_B):
    if arePlotsSaved:
        plt.figure()
        plot_sigma(sigma_BN)
        save_plot("sigma_BN")
        plt.figure()
        plot_omega(omega_BN_B)
        save_plot("omega_BN_B")
    else:
        plt.figure()
        plt.subplot(211)
        plt.title('Sc Att: $\sigma_{BN}$')
        plot_sigma(sigma_BN)
        plt.subplot(212)
        plt.title('Sc Rate: $^B{\omega_{BN}}$')
        plot_omega(omega_BN_B)
    return

def plot_filterOut(sigma_err, omega_err, covar):
    if arePlotsSaved:
        plt.figure()
        plot_sigmaCovar(sigma_err, covar[:,0:4])
        save_plot("sigma_BN error and covariance")
        plt.figure()
        plot_omegaCovar(omega_err, np.stack((covar[:,0],covar[:,4],covar[:,5],covar[:,6]), axis=-1))
        save_plot("omega_BN_B error and covariance")
    else:
        plt.figure()
        plt.subplot(211)
        plt.title('Sc Att Estimation: $\sigma_{BN}$')
        plot_sigmaCovar(sigma_err, covar[:,0:4])
        plt.subplot(212)
        plt.title('Sc Rate Estimation: $^B{\omega_{BN}}$')
        plot_omegaCovar(omega_err, np.stack((covar[:,0],covar[:,4],covar[:,5],covar[:,6]), axis=-1))
    return

def plot_filterPostFits(postFits, noise):
    if arePlotsSaved:
        plt.figure()
        plot_postfit(np.stack((postFits[:,0],postFits[:,1]), axis=-1), noise[1,1])
        save_plot("sigma_BN PostFit First Component")
        plt.figure()
        plot_postfit(np.stack((postFits[:,0],postFits[:,2]), axis=-1), noise[1,1])
        save_plot("sigma_BN PostFit Second Component")
        plt.figure()
        plot_postfit(np.stack((postFits[:,0],postFits[:,3]), axis=-1), noise[1,1])
        save_plot("sigma_BN PostFit Second Component")
    else:
        plt.figure()
        plt.subplot(311)
        plt.title('PostFit $\sigma_{BN}$ First Comp')
        plot_postfit(np.stack((postFits[:,0],postFits[:,1]), axis=-1), noise[1,1])
        plt.subplot(312)
        plt.title('PostFit $\sigma_{BN}$ Second Comp')
        plot_postfit(np.stack((postFits[:,0],postFits[:,2]), axis=-1), noise[1,1])
        plt.subplot(313)
        plt.title('PostFit $\sigma_{BN}$ Third Comp')
        plot_postfit(np.stack((postFits[:,0],postFits[:,3]), axis=-1), noise[1,1])
    return



def plot_orbit(r_BN):
    plt.figure()
    plt.xlabel('$R_x$, km')
    plt.ylabel('$R_y$, km')
    plt.plot(r_BN[:, 1]*m2km, r_BN[:, 2]*m2km, color_x)
    plt.scatter(0, 0)
    if arePlotsSaved:
        save_plot("sc_orbit")
    else:
        plt.title('Spacecraft Orbit')
    return

def plot_multi_orbit_0(dict_data_color):
    plt.figure()
    plt.xlabel('$R_x$, km')
    plt.ylabel('$R_y$, km')
    plt.title('Celestial Bodies IC')
    legend = []
    x_list = []
    y_list = []
    for k, v in dict_data_color.items():
        data = v[0]
        color = v[1]
        legend.append(k)
        x_list.append(np.abs(data[0, 1])) # initial (t=0) x coordinate
        y_list.append(np.abs(data[0, 2])) # initial (t=0) y coordinate
        plt.scatter(data[0, 1], data[0, 2], color=color)
        #plot_celestial_orbit(v[0], v[1])
    x_lim = np.max(x_list)
    y_lim = np.max(y_list)
    f = 1.2
    plt.xlim([-x_lim*f, x_lim*f])
    plt.ylim([-y_lim*f, y_lim*f])
    plt.legend(legend, loc='lower right')
    print 'MultiOrbit legend = ', legend
    plt.axis('equal')
    return

def plot_spacecraft_orbit_0(dict_data_color, r_BN):
    plt.figure()
    plt.xlabel('$R_x$, km')
    plt.ylabel('$R_y$, km')
    plt.title('Spacecraft Cel Body Orbit')
    color_orbit = 'b'
    r_PN_orbit = np.zeros(3)
    r_min = 1E80
    for k, v in dict_data_color.items(): # {planet: [r_PN, 'color']}
        r_PN = v[0]
        r_BP_norm = la.norm(r_BN[0, 1:] - r_PN[0, 1:])
        print k, ': r_BP_norm = ', r_BP_norm
        if r_BP_norm < r_min:
            r_min = r_BP_norm
            r_PN_orbit = r_PN[0, 1:].copy()
            color_orbit = v[1]
    plt.plot(r_BN[:, 1]*m2km, r_BN[:, 2]*m2km, 'magenta')
    plt.scatter(r_PN_orbit[0]*m2km, r_PN_orbit[1]*m2km, color=color_orbit)
    plt.legend(['spacecraft', 'closest planet'])
    plt.axis('equal')
    return


def plot_spacecraft_orbit(dict_data_color, r_BN):
    plt.figure()
    plt.xlabel('$R_x$, km')
    plt.ylabel('$R_y$, km')
    plt.title('Sc and Cel Body Orbit Evolution')
    legend = ['spacecraft']
    plt.plot(r_BN[:, 1]*m2km, r_BN[:, 2]*m2km, 'magenta')
    for k, v in dict_data_color.items():
        r_PN = v[0]
        color = v[1]
        legend.append(k)
        plt.plot(r_PN[:, 1]*m2km, r_PN[:, 2]*m2km, color =color, marker = 'o', linestyle='',linewidth=4)
    plt.legend(legend, loc='lower right')
    plt.axis('equal')








