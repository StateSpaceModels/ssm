{% extends "ordered.tpl" %}

{% block code %}

{% for x in observed %}

static double f_likelihood_tpl_{{ x.name }}(double y, ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double like;
    double *X = p_X->proj;

    {% if x.distribution == 'discretized_normal' %}
    
    double gsl_mu = {{ x.mean }};
    double gsl_sd = {{ x.sd }};


    if (y > 0.0) {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd) - gsl_cdf_gaussian_P(y - 0.5 - gsl_mu, gsl_sd);
    } else {
        like = gsl_cdf_gaussian_P(y + 0.5 - gsl_mu, gsl_sd);
    }

    {% elif x.distribution == 'poisson' %}

    double gsl_mu = {{ x.mean }};

    like = gsl_ran_poisson_pdf(rint(y), gsl_mu);

    {% elif x.distribution == 'binomial' %}

    double gsl_p = {{ x.p }};
    double gsl_n = {{ x.n }};

    like = gsl_ran_binomial_pdf(rint(y), gsl_p, gsl_n);

    {% endif %}

    return like;
}


static double f_obs_mean_tpl_{{ x.name }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return {{ x.mean }};
}

static double f_obs_var_tpl_{{ x.name }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return pow({{ x.sd }}, 2);
}


static double f_obs_ran_tpl_{{ x.name }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;

    {% if x.distribution == 'discretized_normal' %}
    
    double gsl_mu = {{ x.mean }};
    double gsl_sd = {{ x.sd }};

    double yobs = gsl_mu + gsl_ran_gaussian(calc->randgsl, gsl_sd);

    yobs = (yobs >0) ? yobs : 0.0;

    {% elif x.distribution == 'poisson' %}

    double gsl_mu = {{ x.mean }};

    double yobs = gsl_ran_poisson(calc->randgsl, gsl_mu);

    {% elif x.distribution == 'binomial' %}

    double gsl_p = {{ x.p }};
    double gsl_n = {{ x.n }};

    like = gsl_ran_binomial(calc->randgsl, gsl_p, gsl_n);

    {% endif %}

    return yobs;
}

{% endfor %}


/**
 * Approximation of the variance of a function of one random variable
 * First order Taylor expansion
 * Var(f(X))= \sum_i \frac{\partial f(E(X_i))}{\partial x_i}Var(X_i)+\sum_{i\neqj}\frac{\partial f(E(X_i))}{\partial x_i}\frac{\partial f(E(X_j))}{\partial x_ij}Cov(X_i,X_j)
 */
 {% for h, x in h_grads.items() %}
 {% for t, y in x.items() %}
 static double f_var_pred_tpl_{{ y.name }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, double t)
 {
    double res = 0;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&p_X->proj[m], m, m);

    {% for grad_i in y.grads %}
    {% set outer_loop = loop %}
    {% for grad_ii in y.grads %}
    res += {{ grad_i.Cterm }}*{{ grad_ii.Cterm }}*gsl_matrix_get(&Ct.matrix, {{ grad_i.ind }}, {{ grad_ii.ind }});
    {% endfor %}
    {% endfor %}
    
    return res;
}
{% endfor %}
{% endfor %}

ssm_observed_t **_ssm_observed_new(int *observed_length)
{
    *observed_length = {{ observed|length }};

    ssm_observed_t **observed = malloc(*observed_length * sizeof (ssm_observed_t *));
    if (observed == NULL) {
        ssm_print_err("Allocation impossible for ssm_observed_t **");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i=0; i< *observed_length; i++){
        observed[i] = malloc(sizeof (ssm_observed_t));
        if (observed[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_observed_t *");
            exit(EXIT_FAILURE);
        }
    }

    {% for x in observed %}
    observed[{{ loop.index0 }}]->name = strdup("{{ x.name }}");
    observed[{{ loop.index0 }}]->offset = {{ loop.index0 }};
    observed[{{ loop.index0 }}]->f_likelihood = &f_likelihood_tpl_{{ x.name }};
    observed[{{ loop.index0 }}]->f_obs_mean = &f_obs_mean_tpl_{{ x.name }};
    observed[{{ loop.index0 }}]->f_obs_var = &f_obs_var_tpl_{{ x.name }};
    observed[{{ loop.index0 }}]->f_obs_ran = &f_obs_ran_tpl_{{ x.name }};
    observed[{{ loop.index0 }}]->f_var_pred = &f_var_pred_tpl_{{ x.name }};
    {% endfor %}

    return observed;
}


{% endblock %}
