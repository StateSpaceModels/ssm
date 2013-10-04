{% extends "ordered.tpl" %}

{% block code %}

/**
 * Checking if initial conditions are valid
 */
ssm_err_code_t ssm_check_IC(ssm_par_t *par, ssm_calc_t *calc)
{

    ssm_err_code_t cum_status = SSM_SUCCESS;

    /* checking remainders */
    {% for rem, def in f_remainders_par.items() %}
    cum_status |= ({{ def }} < 0.0) ? SSM_ERR_IC : cum_status; //{{ rem }} {% endfor %}

    /* checking state variables for populations without remainder */
    {% for pop in ic %}
    {% for x in pop %}
    cum_status |= ({{ x }} < 0.0) ? SSM_ERR_IC : cum_status;{% endfor %}{% endfor %}

    return cum_status;
}

{% endblock %}
