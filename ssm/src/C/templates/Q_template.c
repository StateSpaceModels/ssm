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
 * Diffusion function for the Extended Kalman Filter
 */
{% for noises_off, tpl in Q.items() %}
void ssm_evalQ_{{ noises_off }}(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc){

    double *X = p_X->proj;
    int i, j;


    {% if is_diff  %}
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    {% endif %}

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    {% if tpl.Q_proc or tpl.Q_obs %}
    double term;

    {% if is_diff  %}
    double diffed[states_diff->length];
    {% endif %}

    {% if tpl.sf %}
    double _sf[{{ step.sf|length }}];
    {% endif %}

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

    {% endif %} // tpl.Q_proc or tpl.Q_obs


    /* caches */
    {% for sf in step.sf %}
    _sf[{{ loop.index0 }}] = {{ sf }};{% endfor %}

    /*
      Q_proc contains only term involving state variables.
    */
    {% if tpl.Q_proc %}
    
    {% for x in tpl.Q_proc %}
    i = {{ x.i }};
    j = {{ x.j }};
    term = {{ x.rate|safe }};

    gsl_matrix_set(Q, i, j, term);
    
    {% if x.i != x.j %}
    gsl_matrix_set(Q, j, i, term);
    {% endif %}
    
    {% endfor %}
    
    {% endif %}

    /*
      Q_obs contains only term involving at least one observed
      variable.
    */
    {% if tpl.Q_obs %}

    {% for x in tpl.Q_obs %}

    {% if x.i.is_obs %}
    i = states_inc->p[{{ x.i }}]->offset;
    {% else %}
    i = states_sv->p[{{ x.i }}]->offset;
    {% endif %}

    {% if x.j.is_obs %}
    j = states_inc->p[{{ x.j }}]->offset;
    {% else %}
    j = states_sv->p[{{ x.j }}]->offset;
    {% endif %}

    gsl_matrix_set(Q, i, j, term);
    
    {% if x.i != x.j %}
    gsl_matrix_set(Q, j, i, term);
    {% endif %}

    {% endfor %}

    {% endif %}
    

    {% if is_diff  %}
    ///////////////////////////////////////////
    // diff term (volatility^2 on diagonal) //
    ///////////////////////////////////////////
    if(is_diff){
	{% for eq in diff %}
	i = it->p[{{ loop.index0 }}]->offset;
	gsl_matrix_set(Q, i, i, pow({{eq}},2));
	{% endfor %}
    }
    {% endif %}

}
{% endfor %}


