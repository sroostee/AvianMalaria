#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

An SIB-model modelling the number of Susceptibles, Infected, and infected with Both malaria strains
individuals over time. 
Inspired by van Baalen (1995) et al. and Alizon (2008) (Multiple infections etc)


Assumptions:

	- Within-host dynamics follow Lotka-Volterra competition

"""
###########################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


###########################		start parameters	###############################

n1_0 = 1
n2_0 = 1
r1 = 0.5
r2 = 0.8
K1= 10
K2 = 10
alpha12 = 0.5
alpha21 = 0.5

#time
ntimepoints = 1000
t = np.linspace(0,20, ntimepoints)


###########################		growth functions	###############################

def LotkaVolterraCompetition(n, t, r1, r2, K1, K2, alpha12, alpha21):
	#calculate malaria growth in host
	n1 = n0[0]
	n2 = n0[1]
	dn1dt = r1*n1*(1-((n1+alpha12*n2)/K1))
	dn2dt = r2*n2*(1-((n2+alpha21*n1)/K2))
	return dn1dt, dn2dt

n0 = (n1_0, n2_0)

out = odeint(LotkaVolterraCompetition, n0, t, args =(r1, r2, K1, K2, alpha12, alpha21))

n1, n2 = out.T

###########################		Plot system	over time	################################ 

plt.plot(t, n1, label="Strain 1")
plt.plot(t, n2, label="Strain 2")
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()

###########################		Plot system	isoclines	################################ 
plt.plot([0,K2/alpha21],[K2,0], color = "blue", label = "isocline for N1")
plt.plot([0,max(K1, K2/alpha21)+2], [0,0], color = "blue")
plt.plot([0,K1],[K1/alpha12,0], color = "orange", label = "isocline for N2")
plt.plot([0,0], [0,max(K2, K1/alpha12)+2], color = "orange")
plt.axis([-1, max(K1, K2/alpha21)+1, -1, max(K2, K1/alpha12)+1])
plt.legend(loc="best")
plt.grid()
plt.show()

############################  some playing around with direction fields		############


# X = np.arange(-10, 10, 1)
# Y = np.arange(-10, 10, 1)
# U, V = np.meshgrid(X**2, Y**2)

# fig, ax = plt.subplots()
# q = ax.quiver(X, Y, U, V)
# ax.quiverkey(q, X=0.3, Y=1.1, U=10,
#              label='Quiver key, length = 10', labelpos='E')

# plt.show()

# xmax = 4.0
# xmin = -xmax
# D = 20
# ymax = 4.0
# ymin = -ymax
# x = np.linspace(xmin, xmax, D)
# y = np.linspace(ymin, ymax, D)
# X, Y = np.meshgrid(x, y)
# deg = np.arctan(Y**2 - X)
# QP = plt.quiver(X,Y,np.cos(deg),np.sin(deg))
# plt.show()


###########################		Population model 	####################################

def a_pop(c_a, n1, n2):
	a = c_a*(n1+n2)
	return a 

def beta_pop(c_b, n1, n2, h):
	beta = (c_b * (n1+n2))/(1+h*(n1+n2))
	return beta

