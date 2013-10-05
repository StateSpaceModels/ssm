{% extends "ordered.tpl" %}

{% block code %}

//stepping functions for ODE and SDEs

{% for noises_off, func in step.func.items() %}
{% if noises_off == 'ode'%}
int ssm_step_ode(double t, const double X[], double f[], void *params)
{% else %}
void ssm_step_sde_{{ noises_off }}(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{% endif %}
{

    {% if noises_off == 'ode'%}
    ssm_calc_t *calc = (ssm_calc_t *) params;
    ssm_nav_t *nav = calc->_nav;
    ssm_par_t *par = calc->_par;
    {% else %}
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = calc->y_pred;
    ssm_it_states_t *states_sv = nav->states_sv;
    {% endif %}

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;

    double _r[{{ step.caches|length }}];

    {% if step.sf %}
    double _sf[{{ step.sf|length }}];{% endif %}

    {% for noise in func.proc.noises %}
    double {{ noise }};{% endfor %}

    {% if is_diff %}
    int i;
    double diffed[states_diff->length];
    {% if noises_off != 'ode'%}
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    {% endif %}
    {% endif %}

    {% if is_diff %}
    for(i=0; i<states_diff->length; i++){
        ssm_state_t *p = states_diff->p[i];
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

    /* caches */
    {% for sf in step.sf %}
    _sf[{{ loop.index0 }}] = {{ sf }};{% endfor %}

    {% for cache in step.caches %}
    _r[{{ loop.index0 }}] = {{ cache }};{% endfor %}

    /* noises */
    {% for noise in func.proc.noises %}
    {{ noise }} = sqrt(dt)*gsl_ran_ugaussian(calc->randgsl);{% endfor %}

    /*ODE system*/
    {% for eq in func.proc.system %}
    f[{{eq.index}}] {% if noises_off == 'ode'%}={% else %}= X[{{eq.index}}] + {% endif %} {{ eq.eq }};{% endfor %}

    //TODO: drift of the diffusion
    //for(i=0; i<states_diff->length; i++){
    //    f[states_diff->p[i]->offset] = 0.0;
    //}

    /*compute incidence:integral between t and t+1*/
    {% for eq in func.obs %}
    f[states_inc->p[{{ eq.index }}]->offset] {% if noises_off == 'ode'%}={% else %}= X[states_inc->p[{{ eq.index }}]->offset] + {% endif %} {{ eq.eq }};{% endfor %}

    {% if noises_off == 'ode'%}
    return GSL_SUCCESS;
    {% else %}
    //TODO better constraints
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<states_sv->length; i++){
        X[states_sv->p[i]->offset] =  (f[states_sv->p[i]->offset] < 0.0) ? 0.0 : f[states_sv->p[i]->offset];
    }

    for(i=0; i<states_inc->length; i++){
        X[states_inc->p[i]->offset] = (f[states_inc->p[i]->offset] < 0.0) ? 0.0 : f[states_inc->p[i]->offset];
    }
    {% endif %}

}
{% endfor %}

{% endblock %}

