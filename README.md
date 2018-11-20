# RadioLensfit 

**Visibilities simulation and shape model fitting of radio star-forming galaxies for Radio Weak Lensing**

This is a tool for simulating visibilities of observed SF galaxies by using SKA1-MID configuration, and measuring galaxy ellipticities in the visibility domain adopting  a Bayesian model fitting approach. It uses an exponential Sersic model and works in the visibility domain avoiding Fourier Transform. 

v1.0 - single galaxy at the phase centre (MPI+OpenMP parallelization) 

v2.0 - many galaxies in the field of view: source extraction + fitting of a single galaxy at a time (only OpenMP parallelization)

**Installation**

GSL library is required.

1) Edit the Makefile:
- enable/disable MPI and OpenMP (default: serial)
- enable/disable gridding (default: gridding enabled)
- set the compiler and compilation flags you want to use (default: GNU)

2) Make.

[Usage](https://github.com/marziarivi/RadioLensfit/wiki)

## Citing RadioLensfit
If you use RadioLensfit and find it useful, please consider citing the related paper: 
_Radio Weak Lensing Shear Measurement in the Visibility Domain_

_Part  I._ Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881 - [arXiv:1603.04784](https://arxiv.org/abs/1603.04784)

_Part II._ Rivi M., Miller L., 2018, MNRAS, 476, 2053 - [arXiv:1709.01827](https://arxiv.org/abs/1709.01827)
