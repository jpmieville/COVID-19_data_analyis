#!/usr/bin/env python
"""
If we are interested in exploring the longer-term persistence and endemic dynamics of an infectious disease, then
clearly demographic processes will be important. The simplest and most common way of introducing demography into
the SIR model is to assume there is a natural host “lifespan”, 1/μ years. Then, the rate at which individuals
(in any epidemiological class) suffer natural mortality is given by μ. It is important to emphasize that this
factor is independent of the disease and is not intended to reflect the pathogenicity of the infectious agent.
Historically, it has been assumed that μ also represents the population’s crude birth rate, thus ensuring that
total population size does not change through time (dS/dt + dI/dt + dR/dt = 0). This framework is very much geared
toward the study of human infections in developed nations—our approach would be different if the host population
exhibited its own “interesting” dynamics (as is often the case with wildlife populations; see Chapter 5).

Parameters
μ 	is the per capita death rate, and the population level birth rate.
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals
    together with the probability of transmission.
γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which
    determines the average infectious period.
S(0) 	is the initial proportion of the population that are susceptible.
I(0) 	is the initial proportion of the population that are infectious.
All rates are specified in days.

Requirements.
All parameters must be positive, and S(0)+I(0) ≤ 1

"""

####################################################################
###    This is the PYTHON version of program 2.2 from page 27 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the simple SIR epidemic with equal births and deaths.    #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import numpy as np
import pylab as pl
import scipy.integrate as spi

mu = 1 / (70 * 365.0)
beta = 520 / 365.0
gamma = 1 / 7.0
TS = 1.0
ND = 60 * 365
S0 = 0.1
I0 = 1e-4
R0 = 1 - S0 - I0
INPUT = (S0, I0, R0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((3))
    V = INP
    Y[0] = mu - beta * V[0] * V[1] - mu * V[0]
    Y[1] = beta * V[0] * V[1] - gamma * V[1] - mu * V[1]
    Y[2] = gamma * V[1] - mu * V[2]
    return Y  # For odeint


t_start = 0.0
t_end = ND
t_inc = TS
t_range = np.arange(t_start, t_end + t_inc, t_inc)
RES = spi.odeint(diff_eqs, INPUT, t_range)

print(RES)

# Ploting
pl.subplot(311)
pl.plot(RES[:, 0], "-g", label="Susceptibles")
pl.title("program_2_2.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles")
pl.subplot(312)
pl.plot(RES[:, 1], "-r", label="Infectious")
pl.xlabel("Time")
pl.ylabel("Infectious")
pl.subplot(313)
pl.plot(RES[:, 2], "-k", label="Recovereds")
pl.xlabel("Time")
pl.ylabel("Recovereds")
pl.show()
