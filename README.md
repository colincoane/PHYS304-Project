# PHYS304_Project
Scripts and files used for simulating solar system motion

Scripts written and run in MATLAB

JPL Horizons Ephemeris Data hard coded into functions. Function "solarSystemData.m" contains body initial states and masses, and "solarSystemEndData.m" contains expected final states from horizons results files

Functions with naming scheme "methodIntegrator.m" are used to perform integration

Scripts with naming scheme "methodFinal.m" are used to "methodIntegrator.m" functions as well as perform plotting and analysis of results

To perform simulation, edit save directory for plots and run "methodFinal.m" in a MATLAB environment. Plots will be automatically saved as png files and displayed in the MATLAB window. Position and Velocity errors will be displayed in the command window and can be ccessed in the MATLAB workspace.
