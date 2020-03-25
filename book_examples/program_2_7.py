#!/usr/bin/env python
"""
Although the SIR and SEIR model paradigms are a good approximation to the epidemiological characteristics of many
infectious diseases, such as measles or influenza, other infections have a more complex natural history. As an example
of how such complexities can be accommodated in the model, will we consider infections such as hepatitis B, herpes, or
chickenpox, where a proportion of infected individuals may become chronic carriers, transmitting infection at a low
rate for many years.

For diseases with carrier states, susceptible individuals can be infected by either carriers or acutely infectious
individuals. It is generally assumed that the progress of infection within an individual is independent of their
source of infection; that is, those infected by acutely infectious individuals and those infected by carriers are
indistinguishable. A recently infected individual is acutely (highly) infectious for a given period and then either
recovers completely or moves into the carrier class. Such dynamics lead to the following model:

Equations

Parameters
μ 	is the per capita death rate, and the population level birth rate.
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals
    together with the probability of transmission.
γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which
    determines the average infectious period.
ε	is the proportion reduction in transmission from carriers compared to standard infectious individuals
q	is the proportion of infected individuals that become carriers rather than fully recover
Γ	is the recovery rate associated with carriers; hence the reciprocal (1/Γ) is the average time an individual is in
    the carrier class
S(0) 	is the initial proportion of the population that are susceptible.
I(0)	is the initial proportion of the population that are infectious
C(0)	is the initial proportion of the population that are carriers

All rates are specified in days.

Requirements.
All parameters must be positive, and S(0)+I(0)+C(0) ≤ 1.

"""
####################################################################
###    This is the PYTHON version of program 2.7 from page 44 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SICR which includes a carrier class.		           #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

beta = 0.2
epsilon = 0.1
gamma = 0.01
Gamma = 0.001
mu = 1 / (50 * 365.0)
q = 0.4
S0 = 0.1
I0 = 1e-4
C0 = 1e-3
ND = 60 * 365
TS = 1.0
INPUT = (S0, I0, C0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((3))
    V = INP
    Y[0] = mu - beta * V[0] * (V[1] + epsilon * V[2]) - mu * V[0]
    Y[1] = beta * V[0] * (V[1] + epsilon * V[2]) - gamma * V[1] - mu * V[1]
    Y[2] = q * gamma * V[1] - Gamma * V[2] - mu * V[2]
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
pl.title("Program_2_7.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles")
pl.subplot(312)
pl.plot(RES[:, 1], "-r", label="Infectious")
pl.xlabel("Time")
pl.ylabel("Infected")
pl.subplot(313)
pl.plot(RES[:, 1], "-m", label="Carriers")
pl.xlabel("Time")
pl.ylabel("Carriers")
pl.show()
