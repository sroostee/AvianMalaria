#!/usr/bin/env python3

"""
date: 02/04/2019
author. S.J. Roostee

Reproduction rate of a rare mutant single and double infection model


"""
###########################

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# beta= 0.1 #transmission efficiency variable
# alpha= 0.3 #disease induced mortality variable
# epsilon = np.linspace(0, 5, 100) #host exploitation strategy
# beta_e = beta * epsilon #transmission efficiency (made up relationship)
# alpha_e = alpha * 0.1 * epsilon#disease induced mortality (made up relationship)


t = np.linspace(0,10,100) #time in days
beta_1 = 0.1 #transmission parameter for first infection by first clone
alpha_1 = 0.1 #disease induced mortality for first infection by first clone
mu = 0.1 #host mortality from other causes
beta_2 = 0.2 #transmission parameter of first clone after new infection by other clone
alpha_2 = 0.3 #disease induced mortality of first clone after new infection by other clone
gamma_2 = 0.1 #transmission parameter of second clone after infection
h = 0.1 #the force of infection of the resident 
		#(probabilty per unit time for a host to become infected by the resident strain)

p1_0 = 1 #probability of infection by the first clone at time 0 (infection occurs at t=0)
p2_0 = 0 #probablity of infection by the second clone at time 0

B_f = (beta_1 + h*(beta_2/(mu+alpha_2)))/(mu + alpha_1 + h)	#per host tranmission of first clone
B_s = gamma_2/(mu+alpha_2) #per host transmission of second clone

x = 0.3 #density of susceptible hosts (uninfected)
y = 0.1 #density of singly infected hosts
z = 0 #density of doubly infected hosts

R_0 = B_f*x + B_s*y

def host_infection_prob(y, t, mu, alpha_1, alpha_2, h):
	p1 = y[0]
	p2 = y[1]
	dp1dt = -(mu + alpha_1 + h)*p1
	dp2dt = h*p1 - (mu+alpha_2)*p2
	return dp1dt, dp2dt

y0 = (p1_0, p2_0)

out = odeint(host_infection_prob, y0, t, args =(mu, alpha_1, alpha_2, h))
p1_out, p2_out = out.T

plt.plot(t, p1_out, label="single infection")
plt.plot(t, p2_out, label="double infection")
plt.legend(loc="best")
plt.xlabel("t")
plt.ylabel("probabilty")
plt.grid()
plt.show()

