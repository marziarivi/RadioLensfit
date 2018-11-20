# RadioLensfit 

**Visibilities simulation and shape model fitting of radio star-forming galaxies for Radio Weak Lensing**

This is a MPI+OpenMP code for simulating visibilities of observed SF galaxies by using SKA1-MID configuration, and measuring galaxy ellipticities in the visibility domain adopting  a Bayesian model fitting approach. It uses an exponential Sersic model and works in the visibility domain avoiding Fourier Transform. 

v1.0 - single galaxy at the phase centre (MPI+OpenMP parallelization) 

v1.1 - single galaxy in the field of view, natural gridding

# Installation

GSL library is required.

1) Edit the Makefile:
- enable/disable MPI and OpenMP (default: serial)
- enable/disable gridding (default: gridding enabled)
- set the compiler and compilation flags you want to use (default: GNU)

2) Make.

# Usage

The radio telescope configuration is set at the beginning of the main() function for SKA1-MID. The source catalog is generated according to the galaxy parameters distributions.

<code> RadioLensfit.x \<filename u-coord\> \<filename v-coord\> \<nge\> \<shear1\> \<shear2\> \<flux-cut\> </code>

- nge*10 galaxies will be simulated with flux larger than \<flux-cut\> [in muJy].
- **g** = (g1,g2) is the shear to apply to the ellipticity of each galaxy.
- Coordinate files are assumed to be txt containing one coordinate per line.

The code produces a file, called _ellipticities\<n\>.txt_, for each MPI task (n=0,1,...N) where each row contains the following galaxy data:

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

## Citing RadioLensfit
If you use RadioLensfit and find it useful, please consider citing the related paper: 

_Radio Weak Lensing Shear Measurement in the Visibility Domain I. Methodology_, Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881 - [arXiv:1603.04784](https://arxiv.org/abs/1603.04784)
