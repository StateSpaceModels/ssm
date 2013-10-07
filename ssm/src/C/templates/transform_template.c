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
{% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) and (p.prior.lower != p.prior.upper) %}
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
{% endif %}


{% if 'f_user2par' in p %}
static double f_user2par_tpl_{{ p.id }}(double x, ssm_input_t *par, ssm_calc_t *calc)
{
    return {{ p.f_user2par }};
}
{% endif %}
{% if 'f_par2user' in p %}
static double f_par2user_tpl_{{ p.id }}(double x, ssm_par_t *par, ssm_calc_t *calc)
{
    return {{ p.f_par2user }};
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

ssm_parameter_t **ssm_parameters_new(int *parameters_length)
{
    *parameters_length = {{ pars|length }};

    ssm_parameter_t **parameters;
    parameters = malloc({{ pars|length }} * sizeof (ssm_parameter_t *));
    if (parameters == NULL) {
        ssm_print_err("Allocation impossible for ssm_parameter_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< {{ pars|length }}; i++){
        parameters[i] = malloc(sizeof (ssm_parameter_t));
        if (parameters[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_parameter_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for p in pars %}
    //{{ p.id }}
    parameters[{{ loop.index0 }}]->name = strdup("{{ p.id }}");
    parameters[{{ loop.index0 }}]->offset = {{ loop.index0 }};
    parameters[{{ loop.index0 }}]->offset_theta = -1;

    {% if 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and (p.prior.lower !=0 or p.prior.upper !=1) and (p.prior.lower != p.prior.upper) %}
    parameters[{{ loop.index0 }}]->f = &f_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_inv = &f_inv_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_derivative = &f_der_tpl_{{ p.id }};
    parameters[{{ loop.index0 }}]->f_inv_derivative = &f_inv_der_tpl_{{ p.id }};
    {% elif 'prior' in p and 'lower' in p.prior and p.prior.lower ==0 and 'upper' not in p.prior %}
    parameters[{{ loop.index0 }}]->f = &ssm_f_log;
    parameters[{{ loop.index0 }}]->f_inv = &ssm_f_inv_log;
    parameters[{{ loop.index0 }}]->f_derivative = &ssm_f_der_log;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &ssm_f_inv_der_log;
    {% elif 'prior' in p and 'lower' in p.prior and 'upper' in p.prior and p.prior.lower == 0 and p.prior.upper == 1 %}
    parameters[{{ loop.index0 }}]->f = &ssm_f_logit;
    parameters[{{ loop.index0 }}]->f_inv = &ssm_f_inv_logit;
    parameters[{{ loop.index0 }}]->f_derivative = &ssm_f_der_logit;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &ssm_f_inv_der_logit;
    {% else %}
    parameters[{{ loop.index0 }}]->f = &ssm_f_id;
    parameters[{{ loop.index0 }}]->f_inv = &ssm_f_id;
    parameters[{{ loop.index0 }}]->f_derivative = &ssm_f_id;
    parameters[{{ loop.index0 }}]->f_inv_derivative = &ssm_f_id;
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
    parameters[{{ loop.index0 }}]->f_par2user = &ssm_f_user_par_id;
    parameters[{{ loop.index0 }}]->f_user2par = &ssm_f_user_par_id;
    {% endif %}
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
    states = malloc(({{ states|length }} + {{ sde|length }} + {{ remainders|length }}) * sizeof (ssm_state_t *));
    if (states == NULL) {
        ssm_print_err("Allocation impossible for ssm_state_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< ({{ states|length }} + {{ sde|length }}); i++){
        states[i] = malloc(sizeof (ssm_state_t));
        if (states[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_state_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for p in states %}
    //{{ p }}
    states[{{ loop.index0 }}]->name = strdup("{{ p }}");
    states[{{ loop.index0 }}]->offset = {{ loop.index0 }};

    states[{{ loop.index0 }}]->f = &ssm_f_id;
    states[{{ loop.index0 }}]->f_inv = &ssm_f_id;
    states[{{ loop.index0 }}]->f_derivative = &ssm_f_id;
    states[{{ loop.index0 }}]->f_inv_derivative = &ssm_f_id;

    states[{{ loop.index0 }}]->f_remainder = NULL;

    states[{{ loop.index0 }}]->ic = {% if p in par_sv %}parameters[{{ loop.index0 }}]{% else %}NULL{% endif %};
    {% endfor %}


    {% for p in sde %}
    //{{ p.id }}
    states[{{ loop.index0 + states|length }}]->name = strdup("{{ p.id }}");
    states[{{ loop.index0 + states|length }}]->offset = {{ loop.index0 + states|length }};

    {% if 'transformation' in p %}
    states[{{ loop.index0 + states|length }}]->f = &f_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_inv = &f_inv_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_derivative = &f_der_tpl_skl_{{ p.id }};
    states[{{ loop.index0 + states|length }}]->f_inv_derivative = &f_inv_der_tpl_skl_{{ p.id }};
    {% else %}
    states[{{ loop.index0 + states|length }}]->f = &ssm_f_id;
    states[{{ loop.index0 + states|length }}]->f_inv = &ssm_f_id;
    states[{{ loop.index0 + states|length }}]->f_derivative = &ssm_f_id;
    states[{{ loop.index0 + states|length }}]->f_inv_derivative = &ssm_f_id;
    {% endif %}

    states[{{ loop.index0 + states|length }}]->f_remainder = NULL;

    states[{{ loop.index0 + states|length }}]->ic = {% if 'ic' in p %}parameters[{{ p.ic }}]{% else %}NULL{% endif %};
    {% endfor %}


    {% for p in remainders %}
    //{{ p }}
    states[{{ loop.index0 + states|length + sde|length }}]->name = strdup("{{ p }}");
    states[{{ loop.index0 + states|length + sde|length }}]->offset = {{ loop.index0 }}; //offset restart at 0 for remainders as they are appart in ssm_hat_t and non existent in ssm_X_t.

    states[{{ loop.index0 + states|length + sde|length }}]->f = &ssm_f_id;
    states[{{ loop.index0 + states|length + sde|length }}]->f_inv = &ssm_f_id;
    states[{{ loop.index0 + states|length + sde|length }}]->f_derivative = &ssm_f_id;
    states[{{ loop.index0 + states|length + sde|length }}]->f_inv_derivative = &ssm_f_id;

    states[{{ loop.index0 + states|length + sde|length }}]->f_remainder = &f_remainder_tpl_{{ p }};

    states[{{ loop.index0 + states|length + sde|length }}]->ic = NULL;
    {% endfor %}


    return states;
}


{% endblock %}
