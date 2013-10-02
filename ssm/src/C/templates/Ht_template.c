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



void eval_Ht(ssm_X_t *p_X, ssm_row_t *row, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    double *X = p_X->proj;
    int i, j;
    gsl_matrix *Ht = calc->Ht;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;
    ssm_it_states_t *states_sv = nav->states_sv;

    {% if is_diff  %}
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    {% endif %}

    {% if is_diff  %}
    double diffed[states_diff->length];
    {% endif %}

    //reset Ht.
    gsl_matrix_set_zero(Ht);

    {% if is_diff %}
    for(i=0; i<states_diff->length; i++){
	ssm_parameter_t *p = states_diff->p[i];
	{% if noises_off != 'ode'%}
	if(is_diff){
	    diffed[i] = p->f_inv(X[p->offset]);
	} else {
	    diffed[i] = gsl_vector_get(par, p->ic->offset);
	}
	{% else %}
	diffed[i] = gsl_vector_get(par, p->ic->offset);
	{% endif %}
    }
    {% endif %}

    // Derivatives of observed means against state variables
    {% for Ht_i in Ht.Ht_sv %}
    {% set outer_loop = loop %}
    {% for Ht_ii in Ht_i %}
    gsl_matrix_set(Ht, 
		   states_sv->p[{{ outer_loop.index0 }}]->offset,
		   {{ loop.index0 }},
		   {{ Ht_ii|safe }});
    {% endfor %}
    {% endfor %}

    // Derivatives of observed means against incidence variables
    {% for Ht_i in Ht.Ht_inc %}
    {% set outer_loop = loop %}
    {% for Ht_ii in Ht_i %}
    gsl_matrix_set(Ht, 
		   states_inc->p[{{ outer_loop.index0 }}]->offset,
		   {{ loop.index0 }},
		   {{ Ht_ii|safe }});
    {% endfor %}
    {% endfor %}

    // Derivatives of observed means against diffusing variables
    {% for Ht_i in Ht.Ht_diff %}
    {% set outer_loop = loop %}
    {% for Ht_ii in Ht_i %}
    gsl_matrix_set(Ht, 
		   states_diff->p[{{ outer_loop.index0 }}]->offset,
		   {{ loop.index0 }},
		   {{ Ht_ii|safe }});
    {% endfor %}
    {% endfor %}


    for(i=0; i< m; i++){
	for(j=0; j< row->ts_nonan_length; j++){
	    gsl_matrix_set(Ht,i,j) = gsl_matrix_get(Ht,i,row->observed[j]->index);
	}
    }
}
