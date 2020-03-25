#!/usr/bin/env python
"""
Numerous infectious diseases confer no long-lasting immunity, such as rotaviruses, sexually transmitted infections,
and many bacterial infections. For these diseases, a individuals can be infected multiple times throughout their
lives, with no apparent immunity. Here, we concentrate briefly on this class of models, called SIS because
recovery from infection is followed by an instant return to the susceptible pool.

Equations

By far the most common use of the SIS model  is to capture the dynamics of sexually transmitted infections. Even
without births, this set of equations has an endemic equilibrium as recoverying individuals replenish the pool
of susceptibles.

We note that S+I =1, so in practise the S equation is redundant.

Parameters
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals
together with the probability of transmission.

γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which
    determines the average infectious period.
I(0) 	is the initial proportion of the population that are infectious.

All rates are specified in days.

Requirements.
All parameters must be positive, and I(0) ≤ 1. Note that S=1-I
"""
####################################################################
###    This is the PYTHON version of program 2.5 from page 39 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the simple SIS epidemic without births or deaths.        #
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
I0 = 1e-6
ND = 70
TS = 1.0
INPUT = (1.0 - I0, I0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((2))
    V = INP
    Y[0] = -beta * V[0] * V[1] + gamma * V[1]
    Y[1] = beta * V[0] * V[1] - gamma * V[1]
    return Y  # For odeint


t_start = 0.0
t_end = ND
t_inc = TS
t_range = np.arange(t_start, t_end + t_inc, t_inc)
RES = spi.odeint(diff_eqs, INPUT, t_range)

print(RES)

# Ploting
pl.subplot(211)
pl.plot(RES[:, 0], "-g", label="Susceptibles")
pl.title("Program_2_5.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles")
pl.subplot(212)
pl.plot(RES[:, 1], "-r", label="Infectious")
pl.xlabel("Time")
pl.ylabel("Infectious")
pl.show()
