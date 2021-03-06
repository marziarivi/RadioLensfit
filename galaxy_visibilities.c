/*
 * Copyright (c) 2015 Marzia Rivi
 *
 * This file is part of RadioLensfit.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "galaxy_visibilities.h"
#include "random_gaussian.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// Compute flux independent model galaxy visibilities analitically
// model a galaxy at the phase centre: visibilities are real numbers (smearing? as data vis are shifted at zero...)
void model_galaxy_visibilities_at_zero(unsigned int nchannels, double* spec, double* wavenumbers,
                        double e1, double e2, double scalelength, unsigned long int num_coords,
                        double* uu_metres, double* vv_metres, unsigned long int* count, double* Modvis)
{
    double wavenumber,den,uu,vv,k1,k2,spectra,shape,phase;
    double detA = 1.-e1*e1-e2*e2;
    double scale = scalelength*ARCS2RAD;
    double scale_factor = (scale*scale)/(detA*detA);
    
    double sum = 0.;
    unsigned long int nv = 0;
    
    for (unsigned int ch=0; ch<nchannels; ch++)
    {
        spectra = spec[ch];
        wavenumber = wavenumbers[ch];
        
        for (unsigned long int i = 0; i < num_coords; ++i)
        {
            uu = uu_metres[i];
            vv = vv_metres[i];
            
            k1 = (1.+e1)*uu + e2*vv;
            k2 = e2*uu + (1.-e1)*vv;
            
            den = 1. + scale_factor*wavenumber*wavenumber*(k1*k1+k2*k2);
            shape = spectra/(den*sqrt(den));
            Modvis[nv] = shape;
            
#ifdef GRID
            sum += Modvis[nv]*Modvis[nv]*count[i];
#else
            sum += Modvis[nv]*Modvis[nv];
#endif
            nv++;
        }
    }
    
    // normalise
    sum = sqrt(sum);
    unsigned long int nvis = num_coords*nchannels;
    for (unsigned long int k=0; k<nvis; k++) Modvis[k] /= sum;
}
    
    
// model galaxy at the galaxy position
void model_galaxy_visibilities(unsigned int nchannels, double* spec, double* wavenumbers, double band_factor,
                               double acc_time, double e1, double e2, double scalelength, double l,
                               double m, unsigned long int num_coords, double* uu_metres, double* vv_metres,
                               unsigned long int* count, complexd* Modvis)
{
    double wavenumber,den,uu,vv,k1,k2,spectra,shape,phase,smear;
    double detA = 1.-e1*e1-e2*e2;
    double scale = scalelength*ARCS2RAD;
    double scale_factor = (scale*scale)/(detA*detA);
    
    double sum = 0.;
    unsigned long int nv = 0;
 
    for (unsigned int ch=0; ch<nchannels; ch++)
    {
        spectra = spec[ch];
        wavenumber = wavenumbers[ch];
        
        for (unsigned long int i = 0; i < num_coords; ++i)
        {
          uu = uu_metres[i];
          vv = vv_metres[i];
            
          phase = uu*l+vv*m;
          /*if (phase != 0.)
          {
             smear = band_factor*phase;
             smear = sin(smear)/smear;
          }
          else smear = 1.;
          */
          phase = -wavenumber*phase;
 
          k1 = (1.+e1)*uu + e2*vv;
          k2 = e2*uu + (1.-e1)*vv;
        
          den = 1. + scale_factor*wavenumber*wavenumber*(k1*k1+k2*k2);
          shape = spectra/(den*sqrt(den));
          Modvis[nv].real = shape*cos(phase); //*smear;
          Modvis[nv].imag = shape*sin(phase); //*smear;
 
#ifdef GRID
          sum += (Modvis[nv].real*Modvis[nv].real + Modvis[nv].imag*Modvis[nv].imag)*count[i];
#else
          sum += (Modvis[nv].real*Modvis[nv].real + Modvis[nv].imag*Modvis[nv].imag);
#endif
          nv++;
        }
    }
    
    // normalise
    sum = sqrt(sum);
    unsigned long int nvis = num_coords*nchannels;
    for (unsigned long int k=0; k<nvis; k++) { Modvis[k].real /= sum; Modvis[k].imag /= sum; }
}
    
    
// Compute data galaxy visibilities per channel  
void data_galaxy_visibilities(double spectra, double wavenumber, double band_factor, double acc_time,
                              double e1, double e2, double scalelength, double flux, double l, double m,
                              unsigned long int num_coords, double* uu_metres, double* vv_metres, complexd* vis)
{
        double den,u,v,k1,k2,phase,shape;
        double detA = 1.-e1*e1-e2*e2;
        double scale = scalelength*ARCS2RAD;  // scale in rad
        double scale_factor = (scale*scale)/(detA*detA);
    
        for (unsigned long int i = 0; i < num_coords; ++i)
        {
            u = uu_metres[i];
            v = vv_metres[i];
            
            phase = u*l+v*m;
            /*if (phase !=0.)
            {
               smear = band_factor*phase;
               smear = sin(smear)/smear;
            }
            else smear = 1.;
            */
            phase = -wavenumber*phase;
            
            k1 = (1.+e1)*u + e2*v;
            k2 = e2*u + (1.-e1)*v;
                
            den = 1. + scale_factor*wavenumber*wavenumber*(k1*k1+k2*k2);
            shape = spectra*flux/(den*sqrt(den));
            vis[i].real = shape*cos(phase); //*smear;
            vis[i].imag = shape*sin(phase); //*smear;
        }
        
}


// Add a random Gaussian noise component to the visibilities.
// visibilities array order: ntimes x nbaselines
void add_system_noise(gsl_rng * gen, unsigned int num_baselines, unsigned int num_times, complexd* vis, double* sigma)
{
    double r1, r2, s1, s2, std, vis_freq, mean = 0.0;
    unsigned int i=0,t,b;
    
    for (b = 0; b <num_baselines; b++)
        for (t = 0; t < num_times; t++)
        {
            /* Combine antenna std.devs. to evaluate the baseline std.dev.
             * See Wrobel & Walker (1999) */
            //std = sqrt(s1*s2);
            std = sigma[b];
          
            /* Apply noise */
            r1 = random_gaussian(&r2, gen);
            //i = t*num_baselines + b;
            vis[i].real += r1 * std + mean;
            vis[i].imag += r2 * std + mean;
            i++;
        }
    
}
    

// Shift data visibilities phase to the galaxy position
void data_visibilities_phase_shift(double wavenumber, double l, double m, unsigned long int num_coords,
                                   double* uu_metres, double* vv_metres, complexd* vis)
{
    double u,v,phase;
    
    for (unsigned long int i = 0; i < num_coords; ++i)
    {
        u = uu_metres[i];
        v = vv_metres[i];
        
        phase = u*l+v*m;
        phase *= wavenumber;
 
        vis[i].real *= cos(phase);
        vis[i].imag *= sin(phase);
    }

}
 
    
#ifdef __cplusplus
}
#endif
