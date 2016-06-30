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


//  RadioLensfit.cpp
//
//  argv[1]  name of the file contaning u coordinates
//  argv[2]  name of the file contaning v coordinates
//  argv[3]  number of galaxies
//  argv[4]  applied shear 1st component
//  argv[5]  applied shear 2nd component
//  argv[6]  minimum galaxy flux in muJy


#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <new>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>

#include "datatype.h"
#include "generate_random_values.h"
#include "read_coordinates.h"
#include "random_gaussian.h"
#include "galaxy_visibilities.h"
#include "likelihood.h"
#include "marginalise_r.h"
#include "distributions.h"
#ifdef GRID
#include "evaluate_uv_grid.h"
#endif

using namespace std;

int main(int argc, char *argv[])
{
    int nprocs, rank, num_threads=1;
#ifdef USE_MPI
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_tot = MPI_Wtime();
#else
    nprocs=1;
    rank=0;
    
    clock_t start_tot;
    start_tot = clock();
#endif
#ifdef _OPENMP
#pragma omp parallel
    num_threads = omp_get_num_threads();
#endif
    if (rank==0) cout << "Number of OpenMP threads = " << num_threads << endl;
    
    if (argc < 7)
    {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: RadioLensfit.x <filename u-coord> <filename v-coord> <nge> <shear1> <shear2> <min flux>" << endl;
        cout << "corresponding to ng = nge*10 galaxies, g1 = shear1, g2 = shear2 and minimum galaxy flux [muJy]" << endl;
        exit(EXIT_FAILURE);
    }

    // Initialise input data
    unsigned int num_stations = 197;      // Number of stations
    unsigned int num_channels = 12;        // Number of frequency channels
    unsigned int num_times = 480; //1920; // Number of time samples
    double freq_start_hz = 950e+6;        // Start Frequency, in Hz
    //double freq_inc_hz = 1e+6;          // Frequency increment, in Hz
    double channel_bandwidth_hz = 20e+6; // Frequency channel bandwidth, in Hz
    double ref_frequency_hz = 1.4e+9;  //Reference frequency in Hz at which fluxes are measured
    int time_acc = 60; //15;     // accumulation time
    double efficiency = 0.9;     // system efficiency
    double SEFD_SKA = 400e+6;    // System Equivalent Flux Density (in micro-Jy) of each SKA1 antenna
    double SEFD_MKT = 551e+6; // SEFD of each MeerKat antenna (in micro-Jy)
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;
    if (rank==0)
    {
        cout << "Number baselines: " << num_baselines << endl;
        cout << "Number of time snapshots: " << num_times << endl;
        cout << "Number of channels: " << num_channels << endl;
        cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
        cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
        cout << "Starting frequency (Hz): " << freq_start_hz << endl;
        cout << "Accumulation time (sec): " << time_acc << endl;
    }
    
    double sizeGbytes, totGbytes = 0.;
    double fov = 3600*ARCS2RAD; //1.22*C0/(freq_start_hz*diameter);  // 1 degree field of view in RAD
    printf("field of view: %e [rad] %f [arcsec] \n",fov,fov/(ARCS2RAD));
    
    // Allocate and read uv coordinates ------------------------------------------------------------------------------
    // coordinates in the file are ordered as nbaselines x ntimes
    // coordinates in the array will be ordered as ntimes x nbaselines
    unsigned long int num_coords = num_times * num_baselines;
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    sizeGbytes = 2*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    double len, threshold = 0.; //only uv-points above this threshold [metres] will be used
    // read only the baselines above the threshold and update their number
    double maxB = read_coord_ska(argv[1], argv[2], num_times, &num_baselines, uu_metres, vv_metres, threshold, &len);
    num_coords = num_times * num_baselines;
    
    unsigned long int num_vis  = (unsigned long int) num_channels * num_times * num_baselines;
    if (rank==0)
    {
        cout << "New number baselines: " << num_baselines << endl;
        cout << "Num visibilities: " << num_vis << endl;
    }
    
    // Pre-compute wavenumber and spectral factor for each channel ---------------------------------------------------------------------
    // They corresponds to the central frequency of each channel
    
    double *wavenumbers = new double[num_channels];
    double ch_freq = freq_start_hz + 0.5*channel_bandwidth_hz;
    double *spec = new double[num_channels];
    
    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
        wavenumbers[ch] = 2.0 * PI * ch_freq / C0;
        spec[ch] = pow(ch_freq/ref_frequency_hz,-0.7);
        ch_freq += channel_bandwidth_hz;
    }
 
    
    // Allocate Data Visibilities ------------------------------------------------------------------------------------------
    complexd* visData;
    try
    {
        visData = new complexd[num_vis];
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated original visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    
#ifdef GRID
    // Gridding uv coordinates ----------------------------------------------------------------------------------------
    int grid_size = 800;
    unsigned long int ncells = grid_size*grid_size;
    double* grid_u = new double[ncells];
    double* grid_v = new double[ncells];
    unsigned long int* count = new unsigned long int[ncells];
    sizeGbytes = (2*ncells*sizeof(double)+ncells*sizeof(unsigned long int))/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated grid coordinates and array counter: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    unsigned long int grid_ncoords = evaluate_uv_grid(len, num_coords, uu_metres, vv_metres, grid_size, grid_u, grid_v, count);
    unsigned long int grid_nvis = num_channels*grid_ncoords;
    complexd* grid_visData;
    try
    {
        grid_visData = new complexd[grid_nvis];
        sizeGbytes = grid_nvis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated gridded visibilities: " << grid_nvis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    if (rank==0) cout << "grid length = " << 2*len << ", grid size = " << grid_size << endl;
#endif
    
    // define steps in galaxy scalelength (Ro in ARCSEC) ------------------------------
    double Rmin = 0.3;
    double Rmax = 3.5;
    int numR = 29;
    double* Ro = new double[numR];
    double* rprior = new double[numR];
    Ro[0] = 0.;
    Ro[1] = Rmin;
    
    rprior[0]= 0.;
    int nRo=2;

    while (nRo<numR && Ro[nRo-1] < Rmax)
    {
        // quadratic spacing of samples
        double Rinterval = 0.08 + 0.4*pow( ((Ro[nRo-1]-Rmin)/(Rmax-Rmin)), 2);
        Ro[nRo] = Ro[nRo-1] + Rinterval;
        nRo++;
    }
    
    numR = nRo;
    if (Ro[nRo-1]>Rmax) Rmax=Ro[nRo-1];
    if (rank==0) cout << numR-1 << " samples in galaxy scale-length, " << Rmin << " < r0 < " << Rmax << " arcsec" << endl;

    
    // Allocate Model Visibilities ------------------------------------------------------------------------------------------
    int num_models = numR-1;
    complexd* visMod;
    try
    {
#if defined GRID
        unsigned long int model_ncoords = grid_ncoords;
#else
        unsigned long int model_ncoords = num_coords;
#endif
        visMod = new complexd[num_models*model_ncoords*num_channels];
        sizeGbytes = num_models*model_ncoords*num_channels*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated models: num_models= " << num_models << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }
    if (rank==0) cout << "Total Visibilities GBytes: " << totGbytes << endl;
 

    // Generate fake galaxies ---------------------------------------------------------------------------------------------------
    // shear to be applied
    double g1 = atof(argv[4]);
    double g2 = atof(argv[5]);
    double* sigmab = new double[num_baselines];
    
    // generate intrinsic ellipticity
    int NP = 5;    // 2NP = number of sampled orientations (points on the circle of radius |e|) for each ellipticity module
    int nge = atoi(argv[3]);
    
    //setup random number generator
    const gsl_rng_type * G;
    gsl_rng * gen;
    G = gsl_rng_mt19937;  // Mersenne Twister
    gen = gsl_rng_alloc(G);
    
    unsigned long int seed = random_seed();
    gsl_rng_set(gen,seed);

#ifdef USE_MPI
    int my_gal = nge/nprocs;
    int rem = nge%nprocs;
    if (rem)
        if (rank < rem) my_gal++;
#else
    int my_gal = nge;
#endif
    int mygalaxies = my_gal*2*NP;
    
    // generate ellipticities
    double *ge1 = new double[mygalaxies];
    double *ge2 = new double[mygalaxies];
    generate_ellipticity(gen,my_gal,NP,ge1,ge2);
    
    // generate flux values
    double Fmin = atof(argv[6]);
    double Fmax = 200.;
    const double beta = -1.34; // flux prior: S^beta
    double *gflux = new double[my_gal];
    generate_random_data(gen,my_gal,gflux,Fmin,Fmax,flux_CDF,beta);
    
    // generate scalelength
    double *gscale = new double[my_gal];
    for (unsigned long int g=0; g<my_gal; g++)
    {
        double mu = scale_mean(gflux[g]); //power law relation between flux and scalelength
        generate_random_data(gen,1,&(gscale[g]),Rmin,Rmax,r_CDF,mu);
    }
    
    // Set function to be minimized
    likelihood_params par;
    par.numr = numR;
    par.ro = Ro;
    par.rprior = rprior;
#if defined GRID
    par.ncoords = grid_ncoords;
    par.uu = grid_u;
    par.vv = grid_v;
    par.data = grid_visData;
    par.count = count;
#else
    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.data = visData;
    par.count = 0;
#endif
    par.nchannels = num_channels;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
    par.mod = visMod;
    par.sigma = (SEFD_SKA*SEFD_SKA)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
    if (rank==0) cout << "sigma_vis  = " << sqrt(par.sigma) << " muJy" << endl;
    for (unsigned int b=0; b<num_baselines; b++) sigmab[b] = sqrt(par.sigma);
    
    FILE *pFile;
    char filename[100];
    sprintf(filename,"ellipticities%d.txt",rank);
    pFile = fopen(filename,"w");
    fprintf(pFile, "flux | scale | e1 | m_e1 | err1 | e2 | m_e2 | err2 | 1D var | SNR |   l  |  m  | \n");

    int bad = 0;
    int gal = 0;
    int np_max = 30;  // min number of sampling points with likelihood above 5%ML
    
    gsl_multimin_function minex_post;
    minex_post.n = 3;
    minex_post.f = f_posterior;
    minex_post.params = &par;
    
    gsl_multimin_function minex_func;
    minex_func.n = 2;
    minex_func.f = f_likelihood;
    minex_func.params = &par;
    
    // use Simplex algorithm of Nelder and Mead provided by the GLS library to minimize -log(likelihood)
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
/*  gsl_multimin_fminimizer *p = 0;
    gsl_vector *sy, *y;
    y = gsl_vector_alloc (3);
    sy = gsl_vector_alloc (3);
    p = gsl_multimin_fminimizer_alloc (T, 3);
    */
    
    gsl_multimin_fminimizer *s = 0;
    gsl_vector *sx, *x;
    x = gsl_vector_alloc (2);
    sx = gsl_vector_alloc (2);
    s = gsl_multimin_fminimizer_alloc (T, 2);
    
#ifdef USE_MPI
    double data_time = 0.;
    double fitting_time = 0.;
    double start_data,end_data,start_fitting,end_fitting;
#else
    clock_t data_time = 0;
    clock_t fitting_time = 0;
    clock_t start_data,end_data,start_fitting,end_fitting;
#endif
    
    double l0,m0,flux,scalelength,ee1,ee2,den,SNR_vis;
    double radius,orient;
    complexd z1,z2;
    
    for (unsigned long int g=0; g<my_gal; g++)
    {
        // Data simulation --------------------------------------------------------------------------------------------------------------------------
        // positions in RAD
        radius = gsl_rng_uniform(gen)*0.5*fov;
        orient = gsl_rng_uniform(gen)*2*PI;
        
        l0 = radius*cos(orient);
        m0 = radius*sin(orient);
        
        par.l0 = l0;
        par.m0 = m0;
        
        // flux and scalelength are the same for all galaxies within the ellipticity circle
        flux = gflux[g];
        scalelength = gscale[g];
 
        // set log(prior) for scalelength
        double mu = scale_mean(flux);
        for (int nRo=1; nRo<numR; nRo++)
              rprior[nRo] = rfunc(mu,scale_std,Ro[nRo]);
        
        for (int ang=0; ang<2*NP; ang++)
        {
#ifdef USE_MPI
          start_data = MPI_Wtime();
#else
          start_data = clock();
#endif
          ee1 = ge1[gal];
          ee2 = ge2[gal];
          gal++;

          // apply shear g
          z1.real = ee1+g1;            z1.imag = ee2+g2;         // z1 = e+g
          z2.real = 1.+ee1*g1+ee2*g2;  z2.imag = ee2*g1-ee1*g2;  // z2 = 1+conj(g)*e
          den = z2.real*z2.real+z2.imag*z2.imag;
          ee1 = (z1.real*z2.real + z1.imag*z2.imag)/den;
          ee2 = (z1.imag*z2.real - z1.real*z2.imag)/den;        // e = z1/z2
            
          SNR_vis = 0.;
          
          // generate galaxy visibilities
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (unsigned int ch = 0; ch < num_channels; ch++)
          {
             //double freq_start = freq_start_hz+ch*channel_bandwidth_hz;
             unsigned long int ch_vis = ch*num_coords;
             data_galaxy_visibilities(spec[ch], wavenumbers[ch], par.band_factor, time_acc, ee1, ee2, scalelength,
                                      flux, l0, m0, num_coords, uu_metres, vv_metres, &(visData[ch_vis]));
 
             double SNR_ch = 0.;
             for (unsigned long int vs = ch_vis; vs < ch_vis+num_coords; vs++)
                   SNR_ch += visData[vs].real*visData[vs].real+visData[vs].imag*visData[vs].imag;
             SNR_ch /= par.sigma;
           
             // Add a random Gaussian noise component to the visibilities.
#ifdef _OPENMP
#pragma omp critical
#endif
             {
               SNR_vis += SNR_ch;
               add_system_noise(gen, num_baselines, num_times, &(visData[ch_vis]), sigmab);
             }
#ifdef GRID
             // Phase shift data visibilities (to be done after gridding because real data will be gridded)
             data_visibilities_phase_shift(wavenumbers[ch], l0, m0, num_coords, uu_metres, vv_metres, &(visData[ch_vis]));
              
             // gridding visibilities ----------------------------------------------------------------------------------------------------------------
             unsigned int ch_visgrid = ch*grid_ncoords;
             gridding_visibilities(num_coords,uu_metres,vv_metres,&(visData[ch_vis]),len,grid_size,&(grid_visData[ch_visgrid]),count);
#endif
          }

#ifdef USE_MPI
          end_data = MPI_Wtime();
          data_time += end_data - start_data;
          start_fitting = MPI_Wtime();
#else
          end_data = clock();
          data_time += (end_data - start_data)/CLOCKS_PER_SEC;
          start_fitting = clock();
#endif
            
          // Model fitting ----------------------------------------------------------------------------------------------------------------------------
            
          // Search for the maximum posterior to find starting ellipticity points
          int iter = 0;
          int status;
          double size;
  /*
          // Starting point
          gsl_vector_set (y, 0, 0.01);
          gsl_vector_set (y, 1, 0.01);
          gsl_vector_set (y, 2, exp(mu));
            
          // Set initial step sizes to 0.1
          gsl_vector_set (sy, 0, 0.1);
          gsl_vector_set (sy, 1, 0.1);
          gsl_vector_set (sy, 2, 0.02);
            
          gsl_multimin_fminimizer_set (p, &minex_post, y, sy);
            
          do
          {
                iter++;
                status = gsl_multimin_fminimizer_iterate(p);
                if (status) break;
                
                size = gsl_multimin_fminimizer_size(p);
                status = gsl_multimin_test_size (size, 1e-3);
                
          }
          while (status == GSL_CONTINUE && iter < 50 && p->fval < 0.);
  */
          double start_e1 = 0.; //gsl_vector_get(p->x, 0);
          double start_e2 = 0.; //gsl_vector_get(p->x, 1);
        
          // Search for the maximum likelihood
          gsl_vector_set (x, 0, start_e1);
          gsl_vector_set (x, 1, start_e2);
          gsl_vector_set_all (sx, 0.1);
        
          gsl_multimin_fminimizer_set (s, &minex_func, x, sx);
          iter = 0;

          do
          {
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);
            
                if (status) break;
            
                size = gsl_multimin_fminimizer_size(s);
                status = gsl_multimin_test_size (size, 1e-3);
           }
           while (status == GSL_CONTINUE && iter < 50 && s->fval < 0.);
        
           double mes_e1, mes_e2, maxL;
           mes_e1 = gsl_vector_get(s->x, 0);
           mes_e2 = gsl_vector_get(s->x, 1);
           maxL= -s->fval;
           cout << "rank:" << rank << " n. " << gal << " flux = " << flux << " scalelength = " << scalelength << " position [arcsec] (" << l0/(ARCS2RAD) << "," << m0/(ARCS2RAD) << "): Maximum log likelihood = " << maxL << " n.iter = " << iter << " for e = " << mes_e1 << "," << mes_e2 <<  "  original e = " << ee1 << "," << ee2 << endl;
            
           // Likelihood sampling to compute mean and variance
           double var_e1, var_e2, cov_e;
           int error = likelihood_sampling(rank,&mes_e1, &mes_e2, maxL, &par, np_max, &var_e1, &var_e2, &cov_e);
           double oneDimvar = sqrt(var_e1*var_e2-cov_e*cov_e);
           fprintf(pFile, "%f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f | %f \n",flux,scalelength,ee1,mes_e1,sqrt(var_e1), ee2,mes_e2,sqrt(var_e2),oneDimvar,sqrt(SNR_vis),l0/(ARCS2RAD),m0/(ARCS2RAD));
           if (error)
           {
              cout << "ERROR likelihood sampling!" << endl;
              bad++;
           }
#ifdef USE_MPI
           end_fitting = MPI_Wtime();
           fitting_time += end_fitting - start_fitting;
#else
           end_fitting = clock();
           fitting_time += (end_fitting - start_fitting)/CLOCKS_PER_SEC;
#endif
       }
    }

    /*
    gsl_vector_free(y);
    gsl_vector_free(sy);
    gsl_multimin_fminimizer_free(p);
    */
    
    gsl_vector_free(x);
    gsl_vector_free(sx);
    gsl_multimin_fminimizer_free(s);
    
    gsl_rng_free(gen);
 
#ifdef USE_MPI
    double end_tot = MPI_Wtime();
#else
    double end_tot = clock();
#endif

    cout << "rank : " << rank << "removed " << bad << " bad data galaxies" << endl << endl;
#ifdef USE_MPI
    double total_time = end_tot - start_tot;
#else
    clock_t total_time = (end_tot - start_tot)/CLOCKS_PER_SEC;
#endif
    cout << "rank: " << rank << " set up time (sec): " << total_time - data_time - fitting_time << endl;
    cout << "rank: " << rank << " data generation time (sec): " << data_time << endl;
    cout << "rank: " << rank << " data fitting computation time (sec): " << fitting_time << endl;
    cout << "rank: " << rank << " Total time (sec): " << total_time << endl;
        
    if (pFile != 0) fclose(pFile);

    
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visMod;
    delete[] visData;
    delete[] Ro;
    delete[] rprior;
    delete[] sigmab;
    delete[] ge1;
    delete[] ge2;
    delete[] gflux;
    delete[] gscale;
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] wavenumbers;
    delete[] spec;
#ifdef GRID
    delete[] grid_u;
    delete[] grid_v;
    delete[] grid_visData;
#endif
#ifdef USE_MPI
      MPI_Finalize() ;
#endif
    return 0;
}
