# RadioLensfit 

Visibilities simulation and shape model fitting of radio galaxies.

This is a MPI+OpenMP code for simulating visibilities of one observed galaxy by using SKA1-MID configuration, and measuring the galaxy shape in the visibility domain adopting  a Bayesian model fitting approach. It uses an exponential Sersic model and works in the visibility domain avoiding Fourier Transform.

version 1.0 - single galaxy at the phase centre 
version 1.1 - single galaxy in the field of view, normal gridding

## Installation

GSL library is required.

1) Edit the Makefile:
- enable/disable MPI and openMP (default: serial)
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
Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881
