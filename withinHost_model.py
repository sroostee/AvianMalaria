#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

Modelling the dynamics of within-host competition of malaria parasites in avian malaria.
The dynamics follow Lotka-Volterra competition of two competitors which has the form of:
dNidt = ri*Ni*(1-((Ni+alphaij*Nj)/Ki))
dNjdt = rj*Nj*(1-((Nj+alphaji*Ni)/Kj))

Assumptions:

	- Within-host dynamics follow Lotka-Volterra competition
	- Only two competitors are present

"""
###########################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###########################		start parameters	###############################

n1_0 = 1
n2_0 = 1
r1 = 0.2
r2 = 0.5
K1= 10
K2 = 10
alpha12 = 1.1
alpha21 = 1.1

#time
ntimepoints = 1000
t = np.linspace(0,100, ntimepoints)

###########################		growth functions	###############################

def LotkaVolterraCompetition(n, t, r1, r2, K1, K2, alpha12, alpha21):
	#calculate malaria growth in host
	N1 = n[0]
	N2 = n[1]
	dN1dt = r1*N1*(1-((N1+alpha12*N2)/K1))
	dN2dt = r2*N2*(1-((N2+alpha21*N1)/K2))
	return dN1dt, dN2dt

n0 = (n1_0, n2_0)

out = odeint(LotkaVolterraCompetition, n0, t, args =(r1, r2, K1, K2, alpha12, alpha21))

n1, n2 = out.T

###########	Effect of within-host dynamics on delta and beta in the population model 	###########

c_delta = 0.8
c_beta = 0.02
h = 0.5

def delta_pop(c_delta, n1, n2):
	#calculation on parasite dependent death rate
	# the more copies a host contains the higher the parasite-caused death rate is
	delta = c_delta*(n1+n2)
	return delta 

def beta_pop(c_b, n, h):
	#caculation on parasote dependent transmission rate
	#the more copies a host contains the higher the transmission rate is, but with a max
	beta = (c_b*n)/(1+h*n)
	return beta

###########################		a and beta over time depending on n1 and n2 	#######

delta_n1 = delta_pop(c_delta, n1, 0)
delta_n2 = delta_pop(c_delta, 0, n2)
delta_all = delta_pop(c_delta, n1, n2)
beta_n1 = beta_pop(c_beta, n1, h)
beta_n2 = beta_pop(c_beta, n2, h)

if __name__ == "__main__":

###########################		Plot within-host system	over time	################################ 

	plt.plot(t, n1, label="Strain 1")
	plt.plot(t, n2, label="Strain 2")
	plt.legend(loc="best")
	plt.xlabel("t")
	plt.ylabel("number of copies")
	plt.grid()
	plt.show()

	###########################		Plot system	isoclines	########################################## 
	plt.plot([0,K2/alpha21],[K2,0], color = "blue", label = "isocline for N1")
	plt.plot([0,max(K1, K2/alpha21)+2], [0,0], color = "blue")
	plt.plot([0,K1],[K1/alpha12,0], color = "orange", label = "isocline for N2")
	plt.plot([0,0], [0,max(K2, K1/alpha12)+2], color = "orange")
	plt.axis([-1, max(K1, K2/alpha21)+1, -1, max(K2, K1/alpha12)+1])
	plt.legend(loc="best")
	plt.grid()
	plt.show()

	plt.plot(t, delta_n1, label = "delta n1")
	plt.plot(t, delta_n2, label = "delta n2")
	plt.plot(t, delta_all, label = "delta all")
	plt.legend(loc="best")
	plt.grid()
	plt.show()

	plt.plot(t, beta_n1, label = "beta n1")
	plt.plot(t, beta_n2, label = "beta n2")
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