# PTFEsim

PFTE simulation in LAMMPS
09/11/15

Collection of scripts to generate a data file of Nafion 
(and other PTFE membranes) with or without electrodes
for further DPD simulation

DPD parametrisation for Nafion, beads A to P
* A: CF2 CF2 CF2 CF2 CF2 CF2
* B: O CF2 C(CF3)F O CF2
* C: CF3 SO3H
* W: 6 H2O
* E: electrodes from carbon black, 54 carbon atoms per bead (optional)
* P: platinum beads (optional)

Default cell size: a few (5-20) DPD sizes (one = 8.14 AA)


## Workflow
1. Create `input.yaml` file with parameters using `default_ptfe_input.sh`
2. Run `gen_ptfe.py` with appropriate options to create LAMMPS data file
3. Run LAMMPS
Data file contains entries in SI units.


## DL_MESO scripts
Scripts for manipulating IO files for [DL_MESO package](http://www.scd.stfc.ac.uk/SCD/40694.aspx) 
for dissipative particle dynamics.

Generate:
* CONTROL file with given values
* initial configuration of particles (CONFIG file)
* interactions and molecules (FIELDS file)

Systems to simulate:
* binary mixture
* diblock copolymer melt of A/B beads, from [Groot, JCP, 1998](http://dx.doi.org/10.1063/1.476300)
* Nafion membrane with electrodes


## Dependencies
* Numpy, f2py
* `sudo pip install docopt` for nice command line reading
* Compile the Fortran module using
 ```
 $ f2py --fcompiler=gnu95 -c f_rdf.f90 -m f_rdf
 ```

