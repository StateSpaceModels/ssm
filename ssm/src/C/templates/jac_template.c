{% extends "ordered.tpl" %}

{% block code %}

void ssm_eval_jac(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    int i, j;
    gsl_matrix *Ft = calc->_Ft;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;
    ssm_it_states_t *states_sv = nav->states_sv;

    {% if is_diff  %}
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    {% endif %}

    {% if is_diff  %}
    double diffed[states_diff->length];
    {% endif %}

    //some terms are always 0: derivative of the ODE (excluding the observed variable) against the observed variable, derivative of the dynamic of the observed variable against the observed variables, derivative of the drift eq.
    gsl_matrix_set_zero(Ft);

    double _r[{{ jac.caches|length }}];

    {% if jac.sf %}
    double _sf[{{ jac.sf|length }}];{% endif %}

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
    {% for sf in jac.sf %}
    _sf[{{ loop.index0 }}] = {{ sf }};{% endfor %}

    {% for cache in jac.caches %}
    _r[{{ loop.index0 }}] = {{ cache }};{% endfor %}


    //first non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the state variable only ( automaticaly generated code )
    {% for jac_i in jac.jac %}
    {% set outer_loop = loop %}
    {% for jac_ii in jac_i %}
    gsl_matrix_set(Ft,
                   states_sv->p[{{ outer_loop.index0 }}]->offset,
                   states_sv->p[{{ loop.index0 }}]->offset,
                   _r[{{ jac_ii }}]);
    {% endfor %}
    {% endfor %}


    //second non null part of the jacobian matrix: derivative of the dynamic of the observed variable against the state variable only ( automaticaly generated code )
    {% for jac_i in jac.jac_obs %}
    {% set outer_loop = loop %}
    {% for jac_ii in jac_i %}
    gsl_matrix_set(Ft,
                   states_inc->p[{{ outer_loop.index0 }}]->offset,
                   states_sv->p[{{ loop.index0 }}]->offset,
                   _r[{{ jac_ii }}]);
    {% endfor %}
    {% endfor %}



    {% if is_diff %}
    if(is_diff){

        //third non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the diff variable (automaticaly generated code)
        {% for jac_i in jac.jac_diff %}
        {% set outer_loop = loop %}
        {% for jac_ii in jac_i %}
        gsl_matrix_set(Ft,
                       states_sv->p[{{ outer_loop.index0 }}]->offset,
                       states_diff->p[{{ loop.index0 }}]->offset,
                       ssm_diff_derivative(_r[{{ jac_ii.value }}], X, states_diff->p[{{ jac_ii.order }}]));
        {% endfor %}
        {% endfor %}

        //fourth non null part of the jacobian matrix: derivative of obs variable against diff (automaticaly generated code)
        {% for jac_i in jac.jac_obs_diff %}
        {% set outer_loop = loop %}
        {% for jac_ii in jac_i %}
        gsl_matrix_set(Ft,
                       states_inc->p[{{ outer_loop.index0 }}]->offset,
                       nav->states_diff->p[{{ loop.index0 }}]->offset,
                       ssm_diff_derivative(_r[{{ jac_ii.value }}], X, states_diff->p[{{ jac_ii.order }}]));
        {% endfor %}
        {% endfor %}

    }
    {% endif %}

}


{% endblock %}
