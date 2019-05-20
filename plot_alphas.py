#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

Script written to visualise the effect of the Lotka-Volterra competition parameters on the 
number of double infected hosts. 

"""
###########################

import numpy as np
from scipy.integrate import odeint
import pandas as pd
# import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import seaborn as sns
import pandas as pd
import withinHost_model as inHost
import host_model as host

######################### generate alpha values ############################
alpha_ij = np.arange(0.1,1.5,0.1)
alpha_ji = np.arange(0.1,1.5,0.1)

i_ij = []
a_ij = []
a_ji = []

for a_i in alpha_ij:
	for a_j in alpha_ji:
		a_ij.append(a_i)
		a_ji.append(a_j)
		#calculate copies of each strain in the host
		# number of copies of both strains
		n1 = inHost.K1
		n2 = inHost.K2
		#at equilibrium
		n1_12 = (inHost.K1 - a_i * inHost.K2)/(1- a_i*a_j) 
		n2_12 = (inHost.K2 - a_i * inHost.K1)/(1- a_i*a_j) 
		#calculate the number of hosts with a double infection
		host_out = odeint(host.eq_sys, host.y0, host.time, args =(inHost.c_delta1, inHost.c_delta2, 
			inHost.c_beta, inHost.h, n1, n2, n1_12, n2_12, host.mu, host.labda))
		S, I_1, I_2, I_12 = host_out.T
		I_12 = I_12[-1]
		i_ij.append(I_12)


I_alpha_df = pd.DataFrame(dict(alphaij = a_ij, alphaji = a_ji, I12 = i_ij))
I_to_alphas = I_alpha_df.pivot("alphaij", "alphaji", "I12")
print(I_to_alphas)

########################	PLOT 	###########################################

############## heatmap
ax = sns.heatmap(I_to_alphas)
plt.show()

# ############ 3D surface map

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# #Plot the surface.
# surf = ax.plot_surface(ALPHA_ij, ALPHA_ji, I_ij, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# # Customize the z axis.
# ax.set_zlim(--1, )
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()



































############# go to alpha > 1 loop ############

# alphap_ij = np.arange(1.2,1.5,0.1)
# alphap_ji = np.arange(1.2,1.5,0.1)

# ip_ij = []

# print("In the >1 loop now.")

# for ap_i in alphap_ij:
# 	for ap_j in alphap_ji:
# 		withinHost_out = odeint(inHost.LotkaVolterraCompetition, inHost.n0, inHost.t, 
# 			args =(inHost.r1, inHost.r2, inHost.K1, inHost.K2, ap_i, ap_j))
# 		n1, n2 = withinHost_out.T
# 		n1_eq = n1[-1]
# 		n2_eq = n2[-1]
# 		host_out = odeint(host.eq_sys, host.y0, host.time, args =(inHost.c_delta, inHost.c_beta, inHost.h, 
# 			n1_eq, n2_eq, host.mu, host.labda))
# 		S, I_1, I_2, Ip_12 = host_out.T
# 		Ip_12 = Ip_12[-1]
# 		ip_ij.append(Ip_12)

# print("I values for alpha <1 calculated")
# ALPHAp_ij, ALPHAp_ji = np.meshgrid(alphap_ij, alphap_ji)
# ip_ij = np.array(ip_ij)
# Ip_ij = np.reshape(ip_ij, (4,4))

################### go to next plot ####################

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# #Plot the surface.
# surf = ax.plot_surface(ALPHAp_ij, ALPHAp_ji, Ip_ij, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# # Customize the z axis.
# ax.set_zlim(-1, 30)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()
