#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import csv
import ConfigParser
import os
import time

# Time at start
start_time = time.time()
# Simulate one liquor compartment and N wood compartments
# there are Nc components
# +----------+----------------------------+
# | liquor   | wood                       |
# |          |                            |
# |          | 0        | 1      ..      N|
# +----------+----------------------------+
#            +--------------------------> z (spatial dimension)
#            0         dz                 1
# States are number of moles of each reagent in each compartment
#     x[0]   | x[1]    | x[2] ... x[N-1]
#     x[N]  | x[N+1]  |
#       .    | ..
#     x[(N-1)*Nc]| x[NNc-1]
#
# We simulate a reaction
# 1 A -> 1 B
# r1 = kr1*Ca
# 1 B -> 1 C
# r2 = kr2*Cb
# dNdt = S*r*V


def reader(filename):
    """read csv file"""
    dirc = {}
    with open(filename) as f:
        reader = csv.reader(f)
        # First row is headings
        reader.next()
        names, valuestrings, units, descriptions = zip(*list(reader))
    values = [float(s) for s in valuestrings]
    for name in names:
        dirc[name] = values[(names.index(name))]
    return dirc


def reaction_rates(C, x, T):
    """ Calculate reaction rates for a column of component concentrations
    :param C:
    :return: reaction rates
    """
    
    C[C<0] = 0    
    
    CL, CC, CA, CS = C
    
    Nl, Nw = unflatx(x)
    # Get total moles
    mass_frac = Nw.sum(axis=1)*componentsMM/parameters['wood_mass']

    if mass_frac[2] >= parameters['phase_limit_1']:
        kr2 = 0.02
    elif mass_frac[2] >= parameters['phase_limit_2']:
        kr2 = 0.02
    else:
        kr2 = 0.02

    k1 = 36.2*T**0.5*numpy.exp(-4807.69/T)  
    
    dCLdt = k1*CL
    dCCdt = 2.53*(CA**0.11)*dCLdt
    dCAdt = (-4.78e-3*dCLdt + 1.81e-2*dCCdt)*1/1
    
    return numpy.array([dCLdt,
                        dCCdt,
                        dCAdt])    
    
#    return numpy.array([kr1*(CL)*CA,
#                        kr2*CC*CA])


def flatx(liquor, wood):
    """ Return a "flattened" version of the state variables """
    return numpy.concatenate((liquor[:, None], wood), axis=1).flatten()


def unflatx(x):
    """
    :param x: flattened state variables
    :return: liquor, wood state variables reshaped
    """
    rectangle = x.reshape((Ncomponents, parameters['Ncompartments']+1))
    liquor = rectangle[:, 0]  # First column is liquor
    wood = rectangle[:, 1:]  # all the rest are wood
    return liquor, wood


def concentrations(x):
    """ Return concentrations given number of moles vector
    """
    Nl, Nw = unflatx(x)

    # calculate concentrations
    cl = Nl/parameters['liquor_volume']
    cw = Nw/wood_compartment_volume

    return cl, cw


def temp(t):
    """ Temperature function
    """

    T = parameters['Ti'] + t * 0.1
    return T

# Read configuration file
config = ConfigParser.ConfigParser()
configfile = 'config.cfg'

if os.path.exists(configfile):
    config.read('config.cfg')
else:
    message = ("Cannot find config file {0}. "
               "Try copying sample_config.cfg to {0}.").format(configfile)
    raise EnvironmentError(message)

datadir = os.path.expanduser(config.get('paths', 'datadir'))
parameter_filename = os.path.join(datadir, 'parameters.csv')

# Read parameter file
parameters = reader(parameter_filename)

components = ['Lignin', 'Carbohydrate', 'Alkali', 'Sulfur']
# Molar mass
componentsMM = [1., 1., 1., 1.]
Ncomponents = len(components)
# stoicheometric matrix, reagents negative, products positive
S = numpy.array([[-1, 0, 0, 0],
                 [0, -1, 0, 0],
                 [0,  0,-1, 0]]).T
t_end = 100

K = numpy.array([0., 0., 0.1, 0])  # diffusion constant (mol/(m^2.s))
# FIXME: K and D should be specified in a similar way
def D(T):
    
    D_OH = 3.4e-2*(T**0.5)*numpy.exp(-4870/(8.314*T))
    
    return numpy.array([[0.], [0.], [D_OH], [0.]])



total_volume = parameters['liquor_volume'] + parameters['wood_volume']

dz = 1./parameters['Ncompartments']
wood_compartment_volume = parameters['wood_volume']/parameters['Ncompartments']

# Initial conditions
Nliq0 = numpy.array([0., 0., 1., 1.])

Nwood0 = numpy.zeros((Ncomponents, parameters['Ncompartments']))
# Lignin & Carbo content
Nwood0[0, :] = 0.01
Nwood0[1, :] = 0.01
Nwood0[3, :] = 0.1

x0 = flatx(Nliq0, Nwood0)


def dxdt(x, t):
    # assert numpy.all(x>=0)
    T = temp(t)

    # unpack variables
    cl, cw = concentrations(x)

    # All transfers are calculated in moles/second

    # Diffusion between liquor and first wood compartment
    transfer_rate = K*parameters['A']*(cl - cw[:, 0])

    # Flows for each block are due to diffusion
    #                                       v symmetry boundary
    #       +----+      +----+      +----+ ||
    # from  | 0  | d[0] | 1  | d[1] | 2  | ||
    # liq ->|    |----->|    |----->|    |-||
    #       +----+      +----+      +----+ ||

    # diffusion in wood (Fick's law)
    # The last compartment sees no outgoing diffusion due to symmetry
    # FIXME: This calculates gradients for both dimensions
    _, gradcwz = numpy.gradient(cw, dz)
    diffusion = -parameters['A']*D(T)*gradcwz
    diffusion[:, -1] = 0

    # reaction rates in wood
    r = numpy.apply_along_axis(reaction_rates, 0, cw, x, T)
    # change in moles due to reaction
    reaction = S.dot(r)*wood_compartment_volume

    # mass balance for liquor:
    dNliquordt = -transfer_rate
    # in wood, we change due to diffusion (left and right) and reaction
    dNwooddt = reaction - diffusion + numpy.roll(diffusion, 1)
    # plus the extra flow from liquor
    dNwooddt[:, 0] += transfer_rate

    return flatx(dNliquordt, dNwooddt)


def totalmass(x):
    return sum(x)

t = numpy.linspace(0, t_end)
z = numpy.linspace(0, 1, parameters['Ncompartments'])
zl = numpy.array([-dz*2, 0])  # z-coords of liquor
Nt = len(t)

xs, info = scipy.integrate.odeint(dxdt, x0, t, full_output=True)

# Work out concentrations
# TODO: This is probably inefficient
cl, cw = map(numpy.array, zip(*map(concentrations, xs)))

# Time at end of run
print ('Simulation run time: ', time.time() - start_time, 'sec')
# Concentrations
ax = None
cm = plt.get_cmap('cubehelix')
for i, component in enumerate(components):
    ax = plt.subplot(Ncomponents + 1, 1, i+1, sharex=ax)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.pcolormesh(t, zl, numpy.atleast_2d(cl[:, i]), cmap=cm)
    plt.pcolormesh(t, z, cw[:, i, :].T, cmap=cm)
    plt.ylabel('[{}]'.format(component))

# Check that we aren't creating or destroying mass
plt.subplot(Ncomponents+1, 1, Ncomponents+1, sharex=ax)
plt.plot(t, [totalmass(x) for x in xs])
plt.ylabel('Total moles')
plt.ylim(ymin=0)
plt.subplots_adjust(hspace=0)
plt.show()
