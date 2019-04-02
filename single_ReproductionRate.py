#!/usr/bin/env python3

"""
date: 02/04/2019
author. S.J. Roostee

Reproduction rate of a rare mutant single infection model


"""
###########################

import numpy as np
import matplotlib.pyplot as plt

x = 0.3 #density of susceptible hosts
beta= 0.1 #transmission efficiency variable
alpha= 0.3 #disease induced mortality variable
epsilon = np.linspace(0, 5, 100) #host exploitation strategy
beta_e = beta * epsilon #transmission efficiency (made up relationship)
alpha_e = alpha * 0.1 * epsilon#disease induced mortality (made up relationship)
mu = 0.05 #host mortality from other causes

plt.plot(beta_e, (mu + alpha_e))
plt.xlabel("total mortality rate (mu + alpha(epsilon))")
plt.ylabel("transmission efficiency (beta(epsilon))")
plt.show()

#calculate reproduction ratio

R_0 = (beta_e*x)/(mu + alpha_e)

#per host transmission factor

B_e = beta_e/(mu + alpha_e)

#R_0 = B_e*x
