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

#include <stdio.h>
#include <stdlib.h>
#include "evaluate_uv_grid.h"

#ifdef __cplusplus
extern "C" {
#endif

// Compute grid of u,v coordinates (coordinates are put in the center of the cell).
unsigned long int evaluate_uv_grid(double len, unsigned long int ncoords, double* u, double* v, int sizeg, double* grid_u, double* grid_v, unsigned long int* count)
{
    unsigned int i,j;
    unsigned long int p,c;
    unsigned long int size = sizeg*sizeg;
    
    double* temp_grid_u = (double *) malloc(size*sizeof(double));
    double* temp_grid_v = (double *) malloc(size*sizeof(double));
    
    double inc = 2*len/sizeg;
    p = 0;
    for (j = 0; j < sizeg; ++j)
    {
        for (i = 0; i < sizeg; ++i)
        {
          temp_grid_u[p] = (i+0.5)*inc - len;
          temp_grid_v[p] = (j+0.5)*inc - len;
          count[p] = 0;
          p++;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned long int k = 0; k < ncoords; k++)
    {
        unsigned int pu = (unsigned int) ((u[k] + len) / inc);
        unsigned int pv = (unsigned int) ((v[k] + len) / inc);
        unsigned long int pc = (unsigned long int) pv * sizeg + pu;
#ifdef _OPENMP
#pragma omp critical
#endif
        count[pc]++;
    }
 
    c=0;
    for (p = 0; p < size; p++)
    {
        if (count[p])
        {
            grid_u[c] = temp_grid_u[p];
            grid_v[c] = temp_grid_v[p];
            count[c] = count[p];
            c++;
        }
    }
 
    free(temp_grid_u);
    free(temp_grid_v);
 
    
    return c;
}

    
/*
 Compute visibilities at the pg = (ug,vg) grid points
 by adding all the original visibilities at the (u,v)
 points falling in the box centered in pg
 */
void gridding_visibilities(unsigned long int ncoords, double *u, double *v, complexd *vis, double len, int sizeg, complexd *new_vis, unsigned long int *count)
{
    unsigned long int p,c;
    double  inc = 2*len/sizeg;
    unsigned long int size = sizeg*sizeg;
    
    complexd* temp_grid_vis = (complexd *) calloc(size,sizeof(complexd));
    
    for (unsigned long int k = 0; k < ncoords; k++)
    {
        unsigned int pu = (unsigned int) ((u[k] + len) / inc);
        unsigned int pv = (unsigned int) ((v[k] + len) / inc);
        unsigned long int pc = (unsigned long int) pv * sizeg + pu;
 
        temp_grid_vis[pc].real += vis[k].real;
        temp_grid_vis[pc].imag += vis[k].imag;
    }
 
    c=0;
    for (p = 0; p < size; p++)
    {
        if (temp_grid_vis[p].real)
        {
            new_vis[c].real = temp_grid_vis[p].real/count[c];
            new_vis[c].imag = temp_grid_vis[p].imag/count[c];
            c++;
        }
    }
    free(temp_grid_vis);
}
    
    
#ifdef __cplusplus
}
#endif
