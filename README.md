# AvianMalaria
Date: 19th June 2019

# Project Overview:

	This project describes an SIR-based co-infection model. The host dynamics 
	have been linked to the within-host dynamics through the transmission and infection
	dependent death rate. Within-host dynamics are described with Lotka-Volterra competition
	dynamics. 
	Variables values are defined in the script and linked between scripts. This means a change in a
 	withinHost_model variable value affects the host_model output, but can also affect the plot outputs.

## Scripts in this prject:

### Models:
	withinHost_model.py
	invasionAnalysis.py
	host_model.py
		
### Plotting of varying parameter values and effect on host population:
	plot_alphas.py
	plot_K.py
	cdelta_plot.py


# Software:

Python 3.7.1

### packages: 
 
- numpy (version 1.16.2)
- matplotlib (version 3.0.3)
- scipy (version 1.2.1)
- seaborn (version 0.9.0)
- pandas (version 0.24.2)

## Authors:
Suze Roostee (supervised by Jörgen Ripa)
