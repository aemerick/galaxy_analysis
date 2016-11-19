# galaxy_analysis

Analysis code for analyzing Enzo simulations of isolated, idealized (dwarf) galaxies simulated using my 
fork of the Enzo-dev code (https://bitbucket.org/aemerick/enzo-emerick) following star formation with 
individual star particles, multiple feedback sources (massive star and AGB star stellar winds, supernovae,
photoelectric heating, LW radiation, ionizing radiation, and cosmic rays), and tracing stellar yields from
an arbitrary list of individual elemental species. 

Code is still a work in progress, but will be used as an analysis pipeline to remotely process simulation
dumps to (nearly) all information that would be useful for interpretating results and informing future
analysis. Will distill large data sets to more compact HDF5 data sets one could easily offload a cluster
and play around with locally to generate plots. This includes generating the very large array of both
gas and stellar abundance ratio properties given the ~10-15 possible elemental species one could have 
in a simulation.

The first test of the above will be to readily compare the star formation and chemodynamic properties
of a series of simulations varying the underlying feedback physics in an 1.5 pc resolution simulation
of a dwarf galaxy with initial gass mass of around 10^6 - 10^7 solar masses.
