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
 * if jparameters is NULL, do not plug values from datapackage
 */
ssm_input_t *ssm_input_new(json_t *jparameters, ssm_nav_t *nav)
{
    ssm_input_t *input = gsl_vector_calloc(nav->par_all->length);

    {% for p in pars %}
    {% if 'prior' in p and 'distribution' in p.prior and p.prior.distribution == 'fixed' %}
    //{{ p.name }}
    gsl_vector_set(input, {{ order_parameters[p.name] }}, {{ p.prior.value }});
    {%endif%}
    {% endfor %}

    if(!jparameters){
	return input;
    }

    int i, index;
    ssm_it_parameters_t *it = nav->theta_all;

    json_t *resource = json_object_get(jparameters, "resources");

    for(index=0; index< json_array_size(resource); index++){
        json_t *el = json_array_get(resource, index);

        const char* name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "values") == 0) {

            json_t *jvalues = json_object_get(el, "data");

            for(i=0; i<it->length; i++){
                json_t *jval = json_object_get(jvalues, it->p[i]->name);
                if(json_is_number(jval)) {
                    gsl_vector_set(input, it->p[i]->offset, json_number_value(jval));
                } else {
                    char str[SSM_STR_BUFFSIZE];
                    sprintf(str, "error: parameters.values.%s is not a number\n", it->p[i]->name);
                    ssm_print_err(str);
                    exit(EXIT_FAILURE);
                }
            }
            break;
        }
    }

    return input;
}
