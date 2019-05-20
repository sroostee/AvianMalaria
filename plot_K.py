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

######################### generate K values ############################
K_1 = np.arange(5,16,1)
K_2 = np.arange(5,16,1)

i_ij = []
K_1_ij = []
K_2_ij = []

for k1 in K_1:
	for k2 in K_2:
		#store K values for dataframe
		K_1_ij.append(k1)
		K_2_ij.append(k2)
		#calculate copies of each strain in the host
		# number of copies of both strains
		n1 = k1
		n2 = k2
		#at equilibrium
		n1_12 = (k1 - inHost.alpha12 * k2)/(1- inHost.alpha12*inHost.alpha21) 
		n2_12 = (k2 - inHost.alpha21 * k1)/(1- inHost.alpha12*inHost.alpha21) 
		#calculate the number of hosts with a double infection
		host_out = odeint(host.eq_sys, host.y0, host.time, args =(inHost.c_delta1, inHost.c_delta2, 
			inHost.c_beta, inHost.h, n1, n2, n1_12, n2_12, host.mu, host.labda))
		S, I_1, I_2, I_12 = host_out.T
		I_12 = I_12[-1]
		i_ij.append(I_12)


I_K_df = pd.DataFrame(dict(K1 = K_1_ij, K2 = K_2_ij, I12 = i_ij))
I_to_Ks = I_K_df.pivot("K1", "K2", "I12")

########################	PLOT 	###########################################

############## heatmap
ax = sns.heatmap(data = I_to_Ks)
plt.show()

# # ############ 3D surface map

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# #Plot the surface.
# surf = ax.plot_surface(ALPHA_ij, ALPHA_ji, I_ij, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# # Customize the z axis.
# ax.set_zlim(-1, 35)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()





