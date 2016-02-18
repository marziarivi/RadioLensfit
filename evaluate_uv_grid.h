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

#ifndef _evaluate_uv_grid_h
#define _evaluate_uv_grid_h

#ifdef __cplusplus
extern "C" {
#endif
    
#include "datatype.h"
    
unsigned long int evaluate_uv_grid(double len, unsigned long int ncoords, double* u, double* v, int sizeg, double* grid_u, double* grid_v, unsigned long int* count);
    
void gridding_visibilities(unsigned long int ncoords, double *u, double *v, complexd *vis, double len, int sizeg, complexd *new_vis, unsigned long int* count);

    
#ifdef __cplusplus
}
#endif

#endif
