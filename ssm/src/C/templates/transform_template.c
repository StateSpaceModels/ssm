{% extends "ordered.tpl" %}

{% block code %}

{% for p in drifts %}
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
{% if p.prior.distribution == 'uniform' %}
static double f_prior_tpl_{{ p.id }}(double x)
{
    return gsl_ran_flat_pdf(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}
{% elif p.prior.distribution == 'normal' %}
static double f_prior_tpl_{{ p.id }}(double x)
{
    return gsl_ran_gaussian_pdf( (x - {{ p.prior.mean }}), {{ p.prior.sd }} );
}
{% endif %}
{% endif %}

{# we create custom functions for logit_ab transformation (to enclose a and b). This is the case only for logit_ab #}
{% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) %}
static double f_tpl_{{ p.id }}(double x)
{
    return ssm_f_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_inv_tpl_{{ p.id }}(double x)
{
    return ssm_f_inv_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_der_tpl_{{ p.id }}(double x)
{
    return ssm_f_der_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_inv_der_tpl_{{ p.id }}(double x)
{
    return ssm_f_der_inv_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}

static double f_inv_der2_tpl_{{ p.id }}(double x)
{
    return ssm_f_der_inv_logit_ab(x, {{ p.prior.lower }}, {{ p.prior.upper }});
}
{% endif %}


{% if 'f_user2par' in p %}
static double f_user2par_tpl_{{ p.id }}(double x, ssm_input_t *par, ssm_calc_t *calc)
{
    return {{ p.f_user2par }};
}
{% endif %}
{% if 'f_2prior' in p %}
static double f_2prior_tpl_{{ p.id }}(double x, ssm_hat_t *hat, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = hat->states;
    return {{ p.f_2prior }};
}
{% endif %}
{% endfor %}


{% for rem, def in f_remainders.items() %}
static double f_remainder_tpl_{{ rem }}(ssm_X_t *p_X, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return {{ def }};
}
{% endfor %}

{% for rem, var in f_remainders_var.items() %}
static double f_remainder_var_tpl_{{ rem }}(ssm_X_t *p_X, ssm_calc_t *calc, ssm_nav_t *nav, double t)
{
    double *X = p_X->proj;
    int m = nav->states_sv_inc->length + nav->states_diff->length;
    gsl_matrix_const_view Ct = gsl_matrix_const_view_array(&X[m], m, m);	
    return {{ var }};
}
{% endfor %}

ssm_parameter_t **ssm_parameters_new(int *parameters_length)
{
    *parameters_length = {{ pars|length }};

    ssm_parameter_t **parameters;
    parameters = malloc(*parameters_length * sizeof (ssm_parameter_t *));
    if (parameters == NULL) {
        ssm_print_err("Allocation impossible for ssm_parameter_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< *parameters_length; i++){
        parameters[i] = malloc(sizeof (ssm_parameter_t));
        if (parameters[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_parameter_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for p in pars %}
    //{{ p.id }}
    parameters[{{ order_parameters[p.id] }}]->name = strdup("{{ p.id }}");
    parameters[{{ order_parameters[p.id] }}]->offset = {{ order_parameters[p.id] }};
    parameters[{{ order_parameters[p.id] }}]->offset_theta = -1;

    {% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) %}
    parameters[{{ order_parameters[p.id] }}]->f = &f_tpl_{{ p.id }};
    parameters[{{ order_parameters[p.id] }}]->f_inv = &f_inv_tpl_{{ p.id }};
    parameters[{{ order_parameters[p.id] }}]->f_derivative = &f_der_tpl_{{ p.id }};
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative = &f_inv_der_tpl_{{ p.id }};
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative2 = &f_inv_der2_tpl_{{ p.id }};
    {% elif 'prior' in p and 'lower' in p.prior and p.prior.lower ==0 and 'upper' not in p.prior %}
    parameters[{{ order_parameters[p.id] }}]->f = &ssm_f_log;
    parameters[{{ order_parameters[p.id] }}]->f_inv = &ssm_f_inv_log;
    parameters[{{ order_parameters[p.id] }}]->f_derivative = &ssm_f_der_log;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative = &ssm_f_inv_der_log;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative2 = &ssm_f_inv_der2_log;
    {% elif 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and p.prior.lower == 0 and p.prior.upper == 1 %}
    parameters[{{ order_parameters[p.id] }}]->f = &ssm_f_logit;
    parameters[{{ order_parameters[p.id] }}]->f_inv = &ssm_f_inv_logit;
    parameters[{{ order_parameters[p.id] }}]->f_derivative = &ssm_f_der_logit;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative = &ssm_f_inv_der_logit;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative2 = &ssm_f_inv_der2_logit;
    {% else %}
    parameters[{{ order_parameters[p.id] }}]->f = &ssm_f_id;
    parameters[{{ order_parameters[p.id] }}]->f_inv = &ssm_f_id;
    parameters[{{ order_parameters[p.id] }}]->f_derivative = &ssm_f_der_id;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative = &ssm_f_der_id;
    parameters[{{ order_parameters[p.id] }}]->f_inv_derivative2 = &ssm_f_der_id;
    {% endif %}

    {% if 'prior' in p and 'distribution' in p.prior and p.prior.distribution != 'fixed' %}
    parameters[{{ order_parameters[p.id] }}]->f_prior = &f_prior_tpl_{{ p.id }};
    {# TODO: fixed case #}
    {% else %}
    parameters[{{ order_parameters[p.id] }}]->f_prior = NULL;
    {% endif %}
    
    parameters[{{ order_parameters[p.id] }}]->f_user2par = &{% if 'transformation' in p %}f_user2par_tpl_{{ p.id }}{% else %}ssm_f_user_par_id{% endif %};
    parameters[{{ order_parameters[p.id] }}]->f_2prior = &{% if 'f_2prior' in p %}f_2prior_tpl_{{ p.id }}{% else %}ssm_f_2prior_id{% endif %};    
    {% endfor %}

    return parameters;
}


/**
 * Adapt. Here we do it only for diffusions
 */
ssm_state_t **ssm_states_new(int *states_length, ssm_parameter_t **parameters)
{
    *states_length = ({{ states|length }} + {{ sde|length }}  + {{ remainders|length }});

    ssm_state_t **states;
    states = malloc(*states_length * sizeof (ssm_state_t *));
    if (states == NULL) {
        ssm_print_err("Allocation impossible for ssm_state_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< *states_length; i++){
        states[i] = malloc(sizeof (ssm_state_t));
        if (states[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_state_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for p in states %}
    //{{ p }}
    states[{{ order_states[p] }}]->name = strdup("{{ p }}");
    states[{{ order_states[p] }}]->offset = {{ order_states[p] }};

    states[{{ order_states[p] }}]->f = &ssm_f_id;
    states[{{ order_states[p] }}]->f_inv = &ssm_f_id;
    states[{{ order_states[p] }}]->f_derivative = &ssm_f_id;
    states[{{ order_states[p] }}]->f_inv_derivative = &ssm_f_id;

    states[{{ order_states[p] }}]->f_remainder = NULL;

    states[{{ order_states[p] }}]->ic = {% if p in par_sv %}parameters[{{ order_states[p] }}]{% else %}NULL{% endif %};
    {% endfor %}

    {% for p in sde %}
    //{{ p.id }}
    states[{{ order_states['diff__' + p.id] }}]->name = strdup("{{ p.id }}");
    states[{{ order_states['diff__' + p.id] }}]->offset = {{ order_states['diff__' + p.id] }};

    {% if 'transformation' in p %}
    states[{{ order_states['diff__' + p.id] }}]->f = &f_tpl_skl_{{ p.id }};
    states[{{ order_states['diff__' + p.id] }}]->f_inv = &f_inv_tpl_skl_{{ p.id }};
    states[{{ order_states['diff__' + p.id] }}]->f_derivative = &f_der_tpl_skl_{{ p.id }};
    states[{{ order_states['diff__' + p.id] }}]->f_inv_derivative = &f_inv_der_tpl_skl_{{ p.id }};
    {% else %}
    states[{{ order_states['diff__' + p.id] }}]->f = &ssm_f_id;
    states[{{ order_states['diff__' + p.id] }}]->f_inv = &ssm_f_id;
    states[{{ order_states['diff__' + p.id] }}]->f_derivative = &ssm_f_id;
    states[{{ order_states['diff__' + p.id] }}]->f_inv_derivative = &ssm_f_id;
    {% endif %}

    states[{{ order_states['diff__' + p.id] }}]->f_remainder = NULL;

    states[{{ order_states['diff__' + p.id] }}]->ic = {% if 'offset_ic' in p %}parameters[{{ p.offset_ic }}]{% else %}NULL{% endif %};
    {% endfor %}


    {% for p in remainders %}
    //{{ p }}
    states[{{ order_states[p] }}]->name = strdup("{{ p }}");
    states[{{ order_states[p] }}]->offset = {{ loop.index0 }}; //!!!! offset restart at 0 for remainders as they are appart in ssm_hat_t and non existent in ssm_X_t.

    states[{{ order_states[p] }}]->f = &ssm_f_id;
    states[{{ order_states[p] }}]->f_inv = &ssm_f_id;
    states[{{ order_states[p] }}]->f_derivative = &ssm_f_id;
    states[{{ order_states[p] }}]->f_inv_derivative = &ssm_f_id;

    states[{{ order_states[p] }}]->f_remainder = &f_remainder_tpl_{{ p }};

    states[{{ order_states[p] }}]->ic = NULL;
    {% endfor %}

    return states;
}


{% endblock %}
