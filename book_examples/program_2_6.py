#!/usr/bin/env python
"""
We now introduce a refinement to the SIR model (Program 2.2) which takes into account a latent period. The process of
transmission often occurs due to an initial inoculation with a very small number of pathogen units (e.g., a few
bacterial cells or virions). A period of time then ensues during which the pathogen reproduces rapidly within the host,
relatively unchallenged by the immune system. During this stage, pathogen abundance is too low for active transmission
to other susceptible hosts, and yet the pathogen is present. Hence, the host cannot be categorized as susceptible,
infectious, or recovered; we need to introduce a new category for these individuals who are infected but not yet
infectious. These individuals are referred to as Exposed and are represented by the variable E in SEIR models.


Parameters
μ 	is the per capita death rate, and the population level birth rate.
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals
    together with the probability of transmission.
γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which
    determines the average infectious period.
σ	is the rate at which individuals move from the exposed to the infectious classes. Its reciprocal (1/σ) is the
    average latent (exposed) period.

S(0) 	is the initial proportion of the population that are susceptible.
E(0) 	is the initial proportion of the population that are exposed (infected but not infectious)
I(0)	is the initial proportion of the population that are infectious

All rates are specified in days.

Requirements.
All parameters must be positive, and S(0)+E(0)+I(0) ≤ 1.


"""
####################################################################
###    This is the PYTHON version of program 2.6 from page 41 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SEIR epidemic with equal births and deaths.          #
### Note we no-longer explicitly model the recovered class.	       #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

mu = 1 / (70 * 365.0)
beta = 520 / 365.0
sigma = 1 / 14.0
gamma = 1 / 7.0
ND = 60 * 365.0
TS = 1.0
S0 = 0.1
E0 = 1e-4
I0 = 1e-4
INPUT = (S0, E0, I0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((3))
    V = INP
    Y[0] = mu - beta * V[0] * V[2] - mu * V[0]
    Y[1] = beta * V[0] * V[2] - sigma * V[1] - mu * V[1]
    Y[2] = sigma * V[1] - gamma * V[2] - mu * V[2]
    return Y  # For odeint


t_start = 0.0
t_end = ND
t_inc = TS
t_range = np.arange(t_start, t_end + t_inc, t_inc)
RES = spi.odeint(diff_eqs, INPUT, t_range)

Rec = 1.0 - (RES[:, 0] + RES[:, 1] + RES[:, 2])
print(RES)

# Ploting
pl.subplot(311)
pl.plot(RES[:, 0], "-g", label="Susceptibles")
pl.title("Program_2_6.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles")
pl.subplot(312)
pl.plot(RES[:, 1], "-m", label="Exposed")
pl.plot(RES[:, 2], "-r", label="Infectious")
pl.legend(loc=0)
pl.xlabel("Time")
pl.ylabel("Infected")
pl.subplot(313)
pl.plot(Rec, "-k", label="Recovereds")
pl.xlabel("Time")
pl.ylabel("Recovereds")
pl.show()
