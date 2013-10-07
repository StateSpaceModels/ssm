{% extends "ordered.tpl" %}

{% block code %}

{% for h, x in h_grads.items() %}
{% for t, y in x.items() %}

/**
 * Approximation of the variance of a function of one random variable
 * First order Taylor expansion
 * Var(f(X))= \sum_i \frac{\partial f(E(X_i))}{\partial x_i}Var(X_i)+\sum_{i\neqj}\frac{\partial f(E(X_i))}{\partial x_i}\frac{\partial f(E(X_j))}{\partial x_ij}Cov(X_i,X_j)
 */
static double var_f_pred_{{ y.id }}(ssm_X_t *p_X, ssm_par_t *par, ssm_calc_t *nav, ssm_calc_t *calc, double t)
{
    double res = 0;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&p_X->proj[m], m, m);

    {{ y.grads }}
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

{% endblock %}
