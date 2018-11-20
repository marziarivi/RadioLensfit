# RadioLensfit 

**Visibilities simulation and shape model fitting of radio star-forming galaxies for Radio Weak Lensing**

This is a tool for simulating visibilities of observed SF galaxies by using SKA1-MID configuration, and measuring galaxy ellipticities in the visibility domain adopting  a Bayesian model fitting approach. It uses an exponential Sersic model and works in the visibility domain avoiding Fourier Transform.

_v1.0_ - single galaxy at the phase centre (MPI+OpenMP parallelization) 

_v2.0_ - many galaxies in the field of view: source extraction + fitting of a single galaxy at a time (only OpenMP parallelization)

## Installation

GSL library is required.

1) Edit the Makefile:
- enable/disable MPI and OpenMP (default: serial)
- enable/disable gridding (default: gridding enabled)
- set the compiler and compilation flags you want to use (default: GNU)

2) Make.

## Usage
The radio telescope configuration is set at the beginning of the main() function for SKA1-MID with 12 frequency channels and 60s accumulation time.

<code> RadioLensfit.x \<filename u-coord\> \<filename v-coord\> \<nge\> \<shear1\> \<shear2\> \<flux-cut\> </code>

- nge*10 galaxies will be simulated with flux larger than \<flux-cut\> [in muJy]
- g=(shear1,shear2) is the shear to apply to the ellipticity of each galaxy
- Coordinate files are assumed to be txt containing one coordinate per line.

The code produces a file, called "ellipticities\<n\>.txt", for each MPI task (n=0,1,...N) where each row contains the following galaxy data:
>
- flux [muJy]
- scalelength [arcsec]
- e1 (original)
- m_e1 (measured) 
- err1 (e1 measure error) 
- e2 (original) 
- m_e2 (measured)
- err2 (e2 measure error)
- 1D likelihood variance
- SNR
- l (source coordinates with respect to the phase centre) 
- m

## Citing Radiolensfit
If you use RadioLensfit and find it useful, please consider citing the related paper: 
"Radio Weak Lensing Shear Measurement in the Visibility Domain"

Part  I. Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881 - [arXiv:1603.04784](https://arxiv.org/abs/1603.04784)

Part II. Rivi M., Miller L., 2018, MNRAS, 476, 2053 - [arXiv:1709.01827](https://arxiv.org/abs/1709.01827)
