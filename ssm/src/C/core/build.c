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
#include <jansson.h> //json


void ssm_input_free(ssm_input_t *input)
{
    gsl_vector_free(input);
}

ssm_par_t *ssm_par_new(ssm_nav_t *nav)
{
    return gsl_vector_calloc(nav->par_all->length);
}

void ssm_par_free(ssm_par_t *par)
{
    gsl_vector_free(par);
}

ssm_theta_t *ssm_theta_new(ssm_nav_t *nav)
{
    return gsl_vector_calloc(nav->theta_all->length);
}

void ssm_theta_free(ssm_theta_t *theta)
{
    gsl_vector_free(theta);
}

ssm_var_t *ssm_var_new(ssm_nav_t *nav, json_t *parameters)
{
    gsl_matrix *m = gsl_matrix_calloc(nav->theta_all->length, nav->theta_all->length);


    int i,j, index;
    ssm_it_parameters_t *it = nav->theta_all;
    json_t *resource = json_object_get(parameters, "resource");
    
    for(index=0; index< json_array_size(resource); index++){
	json_t *el = json_array_get(resource, index);

	const char* name = json_string_value(json_object_get(el, "name"));
	if (strcmp(name, "covariance") == 0) {

	    json_t *values = json_object_get(el, "data");

	    for(i=0; i<it->length; i++){
		for(j=0; j<it->length; j++){
		    json_t *cov_i = json_object_get(values, it->p[i]->name);
		    if(cov_i){
			json_t *cov_ij = json_object_get(cov_i, it->p[j]->name);
			if(cov_ij){
			    gsl_matrix_set(m, i, j, json_number_value(cov_ij));
			}
		    }
		}
	    }
	    
	    break;
	}
    }

    return m;
}

void ssm_var_free(ssm_var_t *var){
    gsl_matrix_free(var);
}
