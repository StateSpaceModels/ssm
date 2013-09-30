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
 * Alloc memory for the psr implementation
 */
void ssm_alloc_psr(ssm_calc_t *calc)
{
    unsigned int *tab = ssm_u1_new({{ alloc|length }}); 

    /*automaticaly generated code: dimension of prob and inc*/
    {% for x in alloc %}
    tab[ORDER_{{x.state}}] = {{x.nb_reaction}};
    {% endfor %}

    calc->prob = ssm_d2_var_new({{ alloc|length }}, tab);
    calc->inc = ssm_u2_var_new({{ alloc|length }}, tab);

    free(tab);
}


/**
 * stepping functions for Poisson System with stochastic rates (psr)
 */
void ssm_step_psr(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    double *X = p_X->proj;
    double dt = p_X->dt;

    double sum, one_minus_exp_sum;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;

    /*0-declaration of noise terms (if any)*/
    {% for n in white_noise %}
    double {{ n.id }};{% endfor %}

    double _r[{{ step.caches|length }}];
    {% if step.sf %}
    double _sf[{{ step.sf|length }}];{% endif %}

    {% if is_diff %}
    int i;
    double diffed[states_diff->length];
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    {% endif %}

    {% if is_diff %}
    for(i=0; i<states_diff->length; i++){       
	ssm_parameter_t *p = states_diff->p[i];
	if(is_diff){
	    diffed[i] = p->f_inv(X[p->offset]);
	} else {
	    diffed[i] = gsl_vector_get(par, p->ic->offset);
	}
    }
    {% endif %}

    /*1-generate noise increments (if any) (automaticaly generated code)*/
    {% if white_noise %}
    if(nav->noises_off & SSM_NO_WHITE_NOISE){
	{% for n in white_noise %}
	{{ n.id }} = 1.0;{% endfor %}
    } else {
	{% for n in white_noise %}
	{{ n.id }} = gsl_ran_gamma(calc->randgsl, (dt)/ pow(gsl_vector_get(par, ORDER_{{ n.sd }}), 2), pow(gsl_vector_get(par, ORDER_{{ n.sd }}), 2))/dt;{% endfor %}
    }
    {% endif %}

    /*2-generate process increments (automaticaly generated code)*/
    {% for sf in step.sf %}
    _sf[{{ loop.index0 }}] = {{ sf }};{% endfor %}

    {% for cache in step.caches %}
    _r[{{ loop.index0 }}] = {{ cache }};{% endfor %}

    {{ step.code }}

    /*3-multinomial drawn (automaticaly generated code)*/
    {% for draw in psr_multinomial %}
    ssm_ran_multinomial(calc->randgsl, {{ draw.nb_exit }}, (unsigned int) X[ORDER_{{ draw.state }}], calc->prob[ORDER_{{ draw.state }}], calc->inc[ORDER_{{ draw.state }}]);{% endfor %}

    /*4-update state variables (automaticaly generated code)*/
    //use inc to cache the Poisson draw as thew might be re-used for the incidence computation
    {% for draw in step.poisson %}
    {{ draw }};{% endfor %}

    {{ step.update_code }}

    /*compute incidence:integral between t and t+1 (automaticaly generated code)*/

    {% for eq in step_inc %}
    X[states_inc->p[{{ eq.index }}]->offset] = {{ eq.right_hand_side }};{% endfor %}
}
