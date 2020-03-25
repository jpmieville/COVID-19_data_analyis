#!/usr/bin/env python
"""
Consider a “closed population” without demographics (no births, deaths, or migration). The scenario we have in mind is
a large naive population into which a low level of infectious agent is introduced and where the resulting epidemic
occurs sufficiently quickly that demographic processes are not influential. We also assume homogeneous mixing,
whereby intricacies affecting the pattern of contacts are discarded, yielding βSI as the transmission term. Given the
premise that underlying epidemiological rates are constant, we get the following SIR equations:
SIR equations

Parameters
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals together with the probability of transmission.
γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which determines the average infectious period.
S(0) 	is the initial proportion of the population that are susceptible.
I(0) 	is the initial proportion of the population that are infectious.
All rates are specified in days.

Requirements.
All parameters must be positive, and S(0)+I(0) ≤ 1

"""
####################################################################
###    This is the PYTHON version of program 2.1 from page 19 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the simple SIR epidemic without births or deaths.        #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import numpy as np
import pylab as pl
import scipy.integrate as spi

beta = 1.4247
gamma = 0.14286
TS = 1.0
ND = 70.0
S0 = 1 - 1e-6
I0 = 1e-6
INPUT = (S0, I0, 0.0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((3))
    V = INP
    Y[0] = -beta * V[0] * V[1]
    Y[1] = beta * V[0] * V[1] - gamma * V[1]
    Y[2] = gamma * V[1]
    return Y  # For odeint


t_start = 0.0
t_end = ND
t_inc = TS
t_range = np.arange(t_start, t_end + t_inc, t_inc)
RES = spi.odeint(diff_eqs, INPUT, t_range)

print(RES)

# Plotting
pl.subplot(211)
pl.plot(RES[:, 0], "-g", label="Susceptibles")
pl.plot(RES[:, 2], "-k", label="Recovereds")
pl.legend(loc=0)
pl.title("Program_2_1.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles and Recovereds")
pl.subplot(212)
pl.plot(RES[:, 1], "-r", label="Infectious")
pl.xlabel("Time")
pl.ylabel("Infectious")
pl.show()
