{% extends "ordered.tpl" %}

{% block code %}

/**
 * Diffusion function for the Extended Kalman Filter
 */
{% for noises_off, tpl in Q.items() %}
void ssm_eval_Q_{{ noises_off }}(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    int i, j;
    double term;
    gsl_matrix *Q = calc->_Q;

    {% if tpl.Q_inc %}
    ssm_it_states_t *states_inc = nav->states_inc;
    ssm_it_states_t *states_sv = nav->states_sv;
    {% endif %}


    {% if is_diff  %}
    int is_diff = ! (nav->noises_off & SSM_NO_DIFF);
    ssm_it_states_t *states_diff = nav->states_diff;
    {% endif %}

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    {% if tpl.Q_proc or tpl.Q_inc or tpl.Q_sde %}

    {% if is_diff  %}
    double diffed[states_diff->length];
    {% endif %}

    {% if tpl.sf %}
    double _sf[{{ step.sf|length }}];
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

    {% endif %} // tpl.Q_proc or tpl.Q_inc


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
    term = {{ x.term|safe }};

    gsl_matrix_set(Q, i, j, term);

    {% if x.i != x.j %}
    gsl_matrix_set(Q, j, i, term);
    {% endif %}

    {% endfor %}

    {% endif %}

    /*
      Q_inc contains only term involving at least one incidence
      variable.
    */
    {% if tpl.Q_inc %}

    {% for x in tpl.Q_inc %}

    {% if x.i.is_inc %}
    i = states_inc->p[{{ x.i.ind }}]->offset;
    {% else %}
    i = states_sv->p[{{ x.i.ind }}]->offset;
    {% endif %}

    {% if x.j.is_inc %}
    j = states_inc->p[{{ x.j.ind }}]->offset;
    {% else %}
    j = states_sv->p[{{ x.j.ind }}]->offset;
    {% endif %}

    term = {{ x.term|safe }};

    gsl_matrix_set(Q, i, j, term);

    {% if x.i != x.j %}
    gsl_matrix_set(Q, j, i, term);
    {% endif %}

    {% endfor %}

    {% endif %}


    {% if is_diff  %}
    ///////////////
    // diff term //
    ///////////////
    if(is_diff){
        {% for x in tpl.Q_sde %}
        i = states_diff->p[{{ x.i }}]->offset;
        j = states_diff->p[{{ x.j }}]->offset;

        term = {{ x.term }};

        gsl_matrix_set(Q, i, j, term);
        {% if x.i != x.j %}
        gsl_matrix_set(Q, j, i, term);
        {% endif %}

        {% endfor %}
    }
    {% endif %}

}
{% endfor %}


{% endblock %}
