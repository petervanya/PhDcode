# DPDcoeffs

## 1. Generate dissipative particle dynamics interaction coefficients.

* generate reasonable clusters of 6 water molecules
* add random noise to coordinates using `wiggle_water.py`
* parse coordinates from g09 outfiles of relaxed clusters
* create pairs of clusters and perform SP calculations on them at various distances
  ranging from 4 AA (almost overlap) to ca 8 AA (too far)
* parse the energy curves of each pair


## 2. Set and perform DPD simulations in LAMMPS

Simulate DPD binary mixture or diblock copolymer melt, as researched 
in `Groot and Madden, JCP, 1998`.

Scripts:
* `gen_diblock_copolymer.py`
* `gen_binmixt.py`
* `gen_random_poly.py`
* `order_param.py`

## Dependencies
* Python3: `sudo pip3 install numpy docopt`

Compile `f_rdf.f90` using `f2py3 -c f_rdf.f90 -m f_rdf`
