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

double ssm_mif_cooling(ssm_options *opts, int m);
{
    return pow(opts->a, (double) (m-1));
}



/**
 * if x is a prameter and y an initial condition and n the number of data points:
 *  rescale cov(x,y) by 1/sqrt(n)
 *  rescale cov(x,x) by (1/sqrt(n))*(1/sqrt(n))
 *  do not rescale cov(y,y)
 */
void ssm_mif_scale_var(ssm_var_t *var, ssm_data_t *data, ssm_nav_t *nav)
{
    int i, j;
    int offset_i, offset_j;

    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    ssm_it_parameters_t *fls = nav->theta_icsv_icdiff;       //parameters fitted with fixed lag smoothing (fls)
 
    double inv = 1.0/((double) data->n_obs);
    double sqrt_inv = 1.0/sqrt((double) data->n_obs);
    
    //mif, mif terms: rescale by inv
    for(i=0; i<mif->lengh; i++){
	offset_i = mif->p[i]->offset_theta;
	for(j=0; j<mif->lengh; j++){
	    offset_j = mif->p[j]->offset_theta;
	    gsl_matrix_set(var, offset_i, offset_j, gsl_matrix_get(var, offset_i, offset_j) * inv);
	}
    }

    //mif, fls and fls, mif terms: rescale by sqrt_inv
    for(i=0; i<mif->lengh; i++){
	offset_i = mif->p[i]->offset_theta;
	for(j=0; j<fls->lengh; j++){
	    offset_j = fls->p[j]->offset_theta;
	    gsl_matrix_set(var, offset_i, offset_j, gsl_matrix_get(var, offset_i, offset_j) * sqrt_inv);
	    gsl_matrix_set(var, offset_j, offset_i, gsl_matrix_get(var, offset_j, offset_i) * sqrt_inv);
	}
    }

}
