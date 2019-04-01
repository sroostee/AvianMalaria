#!/usr/bin/env python3
"""
date: 01/04/2019
author. S.J. Roostee

A simple SIR-model modelling the number of Susceptibles, Infected and Recovered individuals
over time. 


"""
###########################


import numpy as np
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Population
N = 1000
# Initial number of infected, recovered, and susceptible
I0 = 1
R0 = 0
S0 = N - I0 - R0
#Contact rates
beta = 0.2
gamma = 0.1
#time
ntimepoints = 1000
t = np.linspace(0,100, ntimepoints)

def deriv(y, t, N, beta, gamma):
	S = y[0]
	I = y[1]
	R = y[2]
	dSdt = -(beta*S*I)/N
	dIdt = (beta*S*I)/N - gamma*I
	dRdt = gamma*I
	return dSdt, dIdt, dRdt


y0 = (S0, I0, R0)

out = odeint(deriv, y0, t, args =(N, beta, gamma))
S,I,R = out.T

plt.plot(t, S, label="Susceptibles")
plt.plot(t, I, label="Infected")
plt.plot(t, R, label="Recovered")
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()
