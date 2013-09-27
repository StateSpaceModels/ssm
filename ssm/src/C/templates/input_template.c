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

ssm_input_t *ssm_input_new(json_t *parameters, ssm_nav_t *nav)
{        
    ssm_input_t *input = gsl_vector_calloc(nav->par_all->length);

    {% for p in pars %}
    {% if 'prior' in p and 'distribution' in p.prior and p.prior.distribution == 'fixed' %}
    //{{ p.id }}
    gsl_vector_set(input, {{ loop.index0 }}, {{ p.prior.value }});
    {%endif%}
    {% endfor %}

    int i, index;
    ssm_it_parameters_t *it = nav->theta_all;

    json_t *resource = json_object_get(parameters, "resource");
    
    for(index=0; index< json_array_size(resource); index++){
	json_t *el = json_array_get(resource, index);

	const char* name = json_string_value(json_object_get(el, "name"));
	if (strcmp(name, "values") == 0) {

	    json_t *values = json_object_get(el, "data");

	    for(i=0; i<it->length; i++){
		gsl_vector_set(input, it->p[i]->offset, json_number_value(json_object_get(values, it->p[i]->name)));
	    }	   
	    break;
	}
    }

    return input;
}
