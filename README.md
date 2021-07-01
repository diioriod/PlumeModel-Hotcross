# PlumeModel-Hotcross
Large eddy simulation of hydrothermal plumes in time varying cross flows

This work is part of a publication submitted to JGR-Oceans

Turbulence properties of a deep-sea hydrothermal plume in a time-variable cross-flow: field and model comparisons for Dante in the Main Endeavour Field

by Ian Adams and Daniela Di Iorio

University of Georgia, Department of Marine Sciences, Athens, GA 30602

A large eddy simulation is used to study a high temperature hydrothermal vent plume in a stratified and tidally modulated crossflow, to identify its turbulence and mixing characteristics. The model parameters and source conditions that are comparable to the vertical velocity and refractive index fluctuations, measured 20 m above the Dante sulfide mound in the Main Endeavour vent Field, are a heat transport of 50 MW over a cross sectional area of 4×4.5 m2. With these model source conditions and output results taken at 20 m above the source with 1 Hz sampling, the shear production of turbulent kinetic energy (TKE), the mean and turbulent transport of TKE, and the buoyancy production/dissipation are quantified showing that shear production dominates. Similarly, thermal variance production, and its mean and turbulent transport, is also quantified showing that the advective term dominates. Because of enhanced entrainment of ambient water into the plume during strong crossflows, all mean and turbulent quantities show tidally modulated values. Assuming steady state, the dissipation rates are evaluated. During strong crossflows, the tilting of the vertical velocity contours and isotherms plays a critical role in the stability of the plume and in creating high shear and thermal gradients on the upstream side of the plume center axis. These dissipation rates are used to quantify the refractive index fluctuations, and given the high thermal dissipation quantities,
is the main contributing factor in acoustic forward scatter.

To run the model:
Start with the shell script Hotcross-SAPELO2.sh to see how to compile the code and run

Output data files consist of:
cfl.out  ! maxu, minu, maxv, minv, maxwm, minw, mtime: max and min velocities with their time to make sure the flows are reasonable

DanteExp1.case ! input variables for the model run, case.inc must include this file name

DanteExp1.cone_values2to25.nc   ! 1 sec output for turbulence analysis

DanteExp1.nc  ! 1 hourly output

DanteExp1snapshot.nc ! sample output taken at the top of the hour

DanteExp1.out ! standard output to a file

Hotcross-SAPELO2  ! the execution file

site_stratification_fromobs.inc ! profiles of temperature and salinity and sigma-t

salt.diagnostics  ! diagnostic file for salinity

temp.diagnostics  ! diagnostic file for temperature 

tracer.diagnostics  ! diagnostic file for a tracer concentration
