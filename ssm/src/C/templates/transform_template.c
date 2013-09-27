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

{% for p in skeletons %}
{% if 'transformation' in p %}
static double f_tpl_skl_{{ p.id }}(double x)
{
    return {{ p.f }};
}

static double f_inv_tpl_skl_{{ p.id }}(double x)
{
    return {{ p.f_inv }};
}

static double f_der_tpl_skl_{{ p.id }}(double x)
{
    return {{ p.f_der }};
}

static double f_inv_der_tpl_skl_{{ p.id }}(double x)
{
    return {{ p.f_der_inv }};
}
{% endif %}
{% endfor %}


{% for p in parameters %}
{% if 'prior' in p %}
{% if p.prior.distribution == 'uniform' and (p.prior.lower != p.prior.upper) %}
static double f_prior_tpl_{{ p.id }}(double x)
{
    return gsl_ran_flat_pdf(x, {{ p.prior.lower }}, {{ p.prior.upper }});
} 
{% elif p.prior.distribution == 'normal' and (p.prior.sd != 0.0) %}
static double f_prior_tpl_{{ p.id }}(double x)
{
    return gsl_ran_gaussian_pdf( (x - {{ p.prior.mu }}), {{ p.prior.sd }} );
}
{% endif %}
{% endif %}

{# we create custom functions for logit_ab transformation (to enclose a and b). This is the case only for logit_ab #}
{% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) and (p.prior.lower != p.prior.upper) %}
static double f_tpl_{{ p.id }}(double x)
{
    return f_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_inv_tpl_{{ p.id }}(double x)
{
    return f_inv_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_der_tpl_{{ p.id }}(double x)
{
    return f_der_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_inv_der_tpl_{{ p.id }}(double x)
{
    return f_der_inv_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});    
}
{% endif %}


{% if 'f_user2par' in p %}
static double f_user2par_tpl_{{ p.id }}(double x, ssm_input_t *par, ssm_calc_t *calc)
{
    return {{ p.f_user2par }};    
}
{% endif %}
{% if 'f_par2user' in p %}
static double f_par2user_tpl_{{ p.id }}(double x, ssm_input_t *par, ssm_calc_t *calc)
{
    return {{ p.f_par2user }};    
}
{% endif %}
{% endfor %}


ssm_parameter_t **ssm_parameters_new(){

    ssm_parameter_t **parameters;
    parameters = malloc({{ pars|length }} * sizeof (ssm_parameter_t *));
    if (parameters == NULL) {
        print_err("Allocation impossible for ssm_parameter_t **");
        exit(EXIT_FAILURE);
    }
    
    int i;
    for(i=0; i< {{ pars|length }}; i++){
	parameters[i] = malloc(sizeof (ssm_parameter_t *));
	if (parameters[i] == NULL) {
	    print_err("Allocation impossible for ssm_parameter_t *");
	    exit(EXIT_FAILURE);
	}
    }

    {% for p in pars %}
    //{{ p.id }}
    parameters[{{ loop.index0 }}]->name = {{ p.id }};
    parameters[{{ loop.index0 }}]->offset = {{ loop.index0 }}; 

    {% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) and (p.prior.lower != p.prior.upper) %}
    parameters[{{ loop.index0 }}]->f = &f_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_inv = &f_inv_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_derivative = &f_der_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_inv_derivative = &f_inv_der_tpl_{{ p.id }};
    {% elif 'prior' in p and 'lower' in p.prior and p.prior.lower ==0 and 'upper' not in p.prior %}
    parameters[{{ loop.index0 }}]->f = &f_log;
    parameters[{{ loop.index0 }}]->f_inv = &f_inv_log;
    parameters[{{ loop.index0 }}]->f_derivative = &f_der_log;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &f_inv_der_log;
    {% elif 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and p.prior.lower == 0 and p.prior.upper == 1 %}
    parameters[{{ loop.index0 }}]->f = &f_logit;
    parameters[{{ loop.index0 }}]->f_inv = &f_inv_logit;
    parameters[{{ loop.index0 }}]->f_derivative = &f_der_logit;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &f_inv_der_logit;
    {% else %}
    parameters[{{ loop.index0 }}]->f = &f_id;
    parameters[{{ loop.index0 }}]->f_inv = &f_id;
    parameters[{{ loop.index0 }}]->f_derivative = &f_id;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &f_id;
    {% endif %}

    {% if 'prior' in p and 'distribution' in p.prior and p.prior.distribution != 'fixed' %}
    parameters[{{ loop.index0 }}]->prior = &f_prior_tpl_{{ p.id }};
    {# TODO: fixed case #}
    {% else %}
    parameters[{{ loop.index0 }}]->prior = NULL;
    {% endif %}

    {% if 'transformation' in p %}
    parameters[{{ loop.index0 }}]->f_par2user = &f_user2par_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_user2par = &f_par2user_tpl_{{ p.id }};
    {% else %}
    parameters[{{ loop.index0 }}]->f_par2user = &f_user_par_id;
    parameters[{{ loop.index0 }}]->f_user2par = &f_user_par_id;
    {% endif %}
    {% endfor %}

    return parameters;
}


/**
 * Adapt. Here we do it only for diffusions
 */
ssm_state_t **ssm_states_new(ssm_parameter_t **parameters){
    
    ssm_state_t **states;
    states = malloc(({{ states|length }} + {{ sde|length }}) * sizeof (ssm_states_t *));
    if (states == NULL) {
        print_err("Allocation impossible for ssm_state_t **");
        exit(EXIT_FAILURE);
    }
    
    int i;
    for(i=0; i< ({{ states|length }} + {{ sde|length }}); i++){
	states[i] = malloc(sizeof (ssm_state_t *));
	if (states[i] == NULL) {
	    print_err("Allocation impossible for ssm_state_t *");
	    exit(EXIT_FAILURE);
	}
    }

    {% for p in states %}
    //{{ p }}
    states[{{ loop.index0 }}]->name = {{ p }};
    states[{{ loop.index0 }}]->offset = {{ loop.index0 }}; 

    states[{{ loop.index0 }}]->f = &f_id;
    states[{{ loop.index0 }}]->f_inv = &f_id;
    states[{{ loop.index0 }}]->f_derivative = &f_id;
    states[{{ loop.index0 }}]->f_inv_derivative = &f_id;

    states[{{ loop.index0 }}]->ic = {% if p in par_sv %}parameters[{{ loop.index0 }}]{% else %}NULL{% endif %};   
    {% endfor %}


    {% for p in sde %}
    //{{ p.id }}
    states[{{ loop.index0 + states|length }}]->name = {{ p.id }};
    states[{{ loop.index0 + states|length }}]->offset = {{ loop.index0 + states|length }}; 

    {% if 'transformation' in p %}
    states[{{ loop.index0 + states|length }}]->f = &f_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_inv = &f_inv_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_derivative = &f_der_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_inv_derivative = &f_inv_der_tpl_skl_{{ p.id }};
    {% else %}
    states[{{ loop.index0 + states|length }}]->f = &f_id;
    states[{{ loop.index0 + states|length }}]->f_inv = &f_id;
    states[{{ loop.index0 + states|length }}]->f_derivative = &f_id;
    states[{{ loop.index0 + states|length }}]->f_inv_derivative = &f_id;
    {% endif %}

    states[{{ loop.index0 + states|length }}]->ic = {% if 'ic' in p %}parameters[{{ p.ic }}]{% else %}NULL{% endif %};  
    {% endfor %}

    return states;    
}
