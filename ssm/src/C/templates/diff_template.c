{% extends "ordered.tpl" %}

{% block code %}

void ssm_compute_diff(ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    {% if diff|length %}

    int i;
    ssm_it_states_t *it = nav->states_diff;

    double dt = p_X->dt;

    double _w[it->length];
    for(i=0; i<it->length; i++){
	_w[i] = gsl_ran_ugaussian(calc->randgsl);
    }

    {% for eq in diff %}
    p_X->proj[it->p[{{ loop.index0 }}]->offset] += sqrt(dt)*({{ eq }});{% endfor %}


    {% endif %}
}


{% endblock %}

