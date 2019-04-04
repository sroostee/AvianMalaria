#!/usr/bin/env python3

"""
date: 04/04/2019
author. S.J. Roostee

Transmission dynamics of malaria between malaria mosquitos and humans, 
based on the article by Gupta et al. (1994)



note to self:

https://institutefordiseasemodeling.github.io/Documentation/general/model-sir.html

"""
###########################

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

x_0 = 0 #proportion of population with immunity
x1_0= 0.1 #proportion of population with infection (infectious)
y_y_0 = 0.1 #proportion of infectious vectors ???
h = 3 #related to the decay of immunity (average duration of 1/h)
s = 2 #time variable; rate of contribution to transmission of 
		#parasites and losing the capacity to produce infective gametocytes (average time is 1/s)
mu = 0.2 #natural death rate of the host
m = 2 #the number of female mosquitos per human host
a = 3 #the rate at which mosquitos bite the host
b = 0.2 #the proportion of infected bites on host that produce an infection
r = 1 #recovery rate from infection
mu_m = 0.5 #average vector mortality rate


def malaria_transmission(y, t, a, b, m, h, mu, r, s, mu_m):
	x = y[0]
	x1 = y[1]
	y_y =y[2]
	dxdt = a*b*m*y_y*(1-x) - (h+mu)*x
	dx1dt = a*b*m*y_y*(1-x) - s*x1
	dydt = a*x1*(1-y_y) - mu_m*y_y
	return dxdt, dx1dt, dydt

y0 = (x_0, x1_0, y_y_0)

t = np.linspace(0,10,100) #time in days

out = odeint(malaria_transmission, y0, t, args =(a, b, m, h, mu, r, s, mu_m))
x_out, x1_out, y_out = out.T

plt.plot(t, x_out, label="with immunity")
plt.plot(t, x1_out, label="with infection")
plt.plot(t, y_out, label="y")
plt.legend(loc="best")
plt.xlabel("t")
plt.ylabel("proportion")
plt.grid()
plt.show()