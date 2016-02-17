# RadioLensfit 

Authors: Marzia Rivi (UCL), Lance Miller (University of Oxford)

Visibilities simulation and Bayesian shape model fitting of radio galaxies.

This is a MPI+OpenMP code for simulating visibilities of one observed galaxy by using SKA1-MID configuration, and measuring the galaxy shape in the visibility domain adopting  a Bayesian model fitting approach. It uses an exponential Sersic model and works in the visibility domain avoiding Fourier Transform.
See paper: M. Rivi, L. Miller, S. Makhathini, F. Abdalla, in prep.

version 1.0 - galaxies at the phase centre, uniform gridding

## Installation

GSL library is required.

1) Edit the Makefile:
- set the compiler and compilation flags you want to use (default GNU)
- enable/disable MPI and openMP

2) Make.

## Usage

<code> RadioLensfit.x \<filename u-coord\> \<filename v-coord\> \<nge\> \<shear1\> \<shear2\> \<flux-cut\> </code>

- nge*10 galaxies will be simulated with flux larger than \<flux-cut\> [in muJy]
- g=(shear1,shear2) is the shear to apply to the ellipticity of each galaxy

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


To compute the cosmic shear use the python script to process the ellipticities files: 

<code> shear.py -nf \<number of files\> </code>
