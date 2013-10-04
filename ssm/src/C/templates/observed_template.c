/**************************************************************************
 *    This file is part of ssm.
 *
 *    ssm is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    ssm is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with ssm.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "ssm.h"

{% for x in observed %}

static double (*f_likelihood_tpl_{{ x.id }}) (double y, ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    double gsl_mu = {{ x.pdf.mean }};
    double gsl_sd = {{ x.pdf.sd }};

    if (y > 0.0) {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd)-gsl_cdf_gaussian_P(y-0.5-gsl_mu, gsl_sd);
    } else {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd);
    }

    return sanitize_likelihood(like);
}

static double (*f_obs_mean_tpl_{{ x.id }}) (ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return {{ x.pdf.mean }};
}

static double (*f_obs_var_tpl_{{ x.id }}) (ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    double *X = p_X->proj;
    return pow({{ x.pdf.sd }}, 2);
}

static double (*f_obs_tpl_{{ x.id }}) (ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *calc, double t)
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
