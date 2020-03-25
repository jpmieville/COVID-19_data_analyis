#!/usr/bin/env python
"""
Numerous infectious diseases are associated with a substantial mortality risk. Examples include malaria, measles,
whooping cough, SARS, and dengue fever, among others. How do we explore the consequences of infection-induced
mortality? Specifically, how do we incorporate a mortality probability into the SIR equations? The obvious
approach would be to add a term such as −mI to the basic equation, where m is a per capita disease-induced mortality
rate for infected individuals. However, this may be tricky to interpret biologically or estimate from data. Instead,
it is preferable to think about the probability, ρ, of an individual in the I class dying from the infection before
either recovering or dying from natural causes.

We initial consider the case of frequency-dependent transmission, where as total population size N decreases, due to
disease-induced mortality, there is no change in the  interaction between hosts. To make the dynamics clearer we
switch to considering the number or density (rather than proportion) of individuals.

Equations

Parameters
ρ	is the mortality probability. The probability than an infected individual dies from the disease before recovering.
μ 	is the per capita death rate from natural causes
ν 	is the population level birth rate, we can equate ν/μ with the concept of a carrying capacity
β 	is the transmission rate and incorporates the encounter rate between susceptible and infectious individuals
    together with the probability of transmission.
γ 	is called the removal or recovery rate, though often we are more interested in its reciprocal (1/γ) which
    determines the average infectious period.
X(0) 	is the initial number or density of susceptible individuals.
Y(0) 	is the initial number or density of infectious individuals.
N(0)	is the initial population size
All rates are specified in days.

Requirements.
All parameters must be positive. Remember, X, Y and ν all refer to numbers; ρ ≤ 1 as it is a probability.

"""

####################################################################
###    This is the PYTHON version of program 2.4 from page 36 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR model with a probability of mortality, and	   #
### unequal births and deaths. This code assumes Frequency-        #
### Dependent Transmission        								   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import numpy as np
import pylab as pl
import scipy.integrate as spi

rho = 0.5
nu = mu = 1 / (70 * 365.0)
beta = 520 / 365.0
gamma = 1 / 7.0
TS = 1.0
ND = 1e5
N0 = 1
X0 = 0.2
Y0 = 1e-4
Z0 = N0 - X0 - Y0
INPUT = (X0, Y0, Z0)


def diff_eqs(INP, t):
    """The main set of equations"""
    Y = np.zeros((3))
    V = INP
    Y[0] = mu - beta * V[0] * V[1] / sum(V) - mu * V[0]
    Y[1] = beta * V[0] * V[1] / sum(V) - (gamma + mu) * V[1] / (1 - rho)
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
pl.title("Program_2_4.py")
pl.xlabel("Time")
pl.ylabel("Susceptibles")
pl.subplot(312)
pl.plot(RES[:, 1], "-r", label="Infectious")
pl.xlabel("Time")
pl.ylabel("Infectious")
pl.subplot(313)
pl.plot(RES[:, 2], "-k", label="Recovereds")
pl.plot(sum((RES[:, 0], RES[:, 1], RES[:, 2])), "--k", label="Total Population")
pl.xlabel("Time")
pl.legend(loc=0)
pl.ylabel("Recovereds\nTotal Population")
pl.show()
