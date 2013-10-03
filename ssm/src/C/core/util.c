/**************************************************************************
 *    This file is part of ssm.
 *
 *    ssm is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    ssm is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with ssm.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/


#include "ssm.h"

/**
 * check if name is in it
 */
int ssm_in_par(ssm_it_parameters_t *it, const char *name)
{
    int i;
    for(i=0; i<it->length; i++){
        if(strcmp(it->p[i]->name, name) == 0){
            return 1;
        }
    }

    return 0;
}



/**
 * tranform --interpolation argument into gsl_interp_type *.
 */
const gsl_interp_type *ssm_str_to_interp_type(const char *optarg){

    if (strcmp(optarg, "linear") == 0) {
        return gsl_interp_linear;
    } else if (strcmp(optarg, "polynomial") == 0){
        return gsl_interp_polynomial;
    } else if (strcmp(optarg, "cspline") == 0){
        return gsl_interp_cspline;
    } else if (strcmp(optarg, "cspline_periodic") == 0){
        return gsl_interp_cspline_periodic;
    } else if (strcmp(optarg, "akima") == 0){
        return gsl_interp_akima;
    } else if (strcmp(optarg, "akima_periodic") == 0){
        return gsl_interp_akima_periodic;
    }

    print_warning("Unknown gsl interpolator for metadata. Linear interpolator will be used instead.");
    return gsl_interp_linear;
}



/**
 * make sure that n_threads <= J and return safe n_threads
 */

int ssm_sanitize_n_threads(int n_threads, ssm_fitness_t *fitness)
{
    if(n_threads > fitness->J){
        return fitness->JJ;
    } else {
        return n_threads;
    }
}


