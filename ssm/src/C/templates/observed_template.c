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

static double f_obs_tpl_{{ x.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    double gsl_mu = {{ x.pdf.mean }};
    double gsl_sd = {{ x.pdf.sd }};

    double yobs= gsl_mu+gsl_ran_gaussian(calc->randgsl, gsl_sd);

    return (yobs >0) ? yobs : 0.0;
}

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
    observed[{{ loop.index0 }}]->f_obs = &f_obs_tpl_{{ x.id }};
    {% endfor %}

    return observed;
}


{% endblock %}
