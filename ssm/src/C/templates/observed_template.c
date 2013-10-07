{% extends "ordered.tpl" %}

{% block code %}

{% for x in observed %}

static double f_likelihood_tpl_{{ x.id }}(double y, ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double like;
    double *X = p_X->proj;
    double gsl_mu = {{ x.pdf.mean }};
    double gsl_sd = {{ x.pdf.sd }};

    if (y > 0.0) {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd) - gsl_cdf_gaussian_P(y - 0.5 - gsl_mu, gsl_sd);
    } else {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd);
    }

    return like;
}

static double f_obs_mean_tpl_{{ x.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return {{ x.pdf.mean }};
}

static double f_obs_var_tpl_{{ x.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return pow({{ x.pdf.sd }}, 2);
}

static double f_obs_ran_tpl_{{ x.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    double gsl_mu = {{ x.pdf.mean }};
    double gsl_sd = {{ x.pdf.sd }};

    double yobs= gsl_mu+gsl_ran_gaussian(calc->randgsl, gsl_sd);

    return (yobs >0) ? yobs : 0.0;
}

{% endfor %}


/**
 * Approximation of the variance of a function of one random variable
 * First order Taylor expansion
 * Var(f(X))= \sum_i \frac{\partial f(E(X_i))}{\partial x_i}Var(X_i)+\sum_{i\neqj}\frac{\partial f(E(X_i))}{\partial x_i}\frac{\partial f(E(X_j))}{\partial x_ij}Cov(X_i,X_j)
 */
{% for h, x in h_grads.items() %}
{% for t, y in x.items() %}
static double var_f_pred_tpl_{{ y.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *nav, ssm_calc_t *calc, double t)
{
    double res = 0;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&p_X->proj[m], m, m);

    {% for grad_i in y.grads %}
    {% set outer_loop = loop %}
    {% for grad_ii in y.grads %}
    res += {{ grad_i.Cterm }}*{{ grad_ii.Cterm }}*gsl_matrix_get(Ct,{{ grad_i.ind }},{{ grad_ii.ind }});
    {% endfor %}
    {% endfor %}
    
    return res;
}
{% endfor %}
{% endfor %}

ssm_observed_t **ssm_observed_new(int *observed_length)
{
    *observed_length = {{ observed|length }};

    ssm_observed_t **observed;
    observed = malloc({{ observed|length }} * sizeof (ssm_observed_t *));
    if (observed == NULL) {
        ssm_print_err("Allocation impossible for ssm_observed_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< {{ observed|length }}; i++){
        observed[i] = malloc(sizeof (ssm_observed_t));
        if (observed[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_observed_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for x in observed %}
    observed[{{ loop.index0 }}]->name = strdup("{{ x.id }}");
    observed[{{ loop.index0 }}]->offset = {{ loop.index0 }};
    observed[{{ loop.index0 }}]->f_likelihood = &f_likelihood_tpl_{{ x.id }};
    observed[{{ loop.index0 }}]->f_obs_mean = &f_obs_mean_tpl_{{ x.id }};
    observed[{{ loop.index0 }}]->f_obs_var = &f_obs_var_tpl_{{ x.id }};
    observed[{{ loop.index0 }}]->f_obs_ran = &f_obs_ran_tpl_{{ x.id }};
    observed[{{ loop.index0 }}]->var_f_pred = &var_f_pred_tpl_{{ x.id }};
    {% endfor %}

    return observed;
}


{% endblock %}
