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

/**
 * simplex algo using GSL. Straightforward adaptation of the GSL doc
 * example
 */
double ssm_simplex(ssm_theta_t *theta, ssm_var_t *var, void *params, double (*f_simplex)(const gsl_vector *x, void *params), ssm_nav_t *nav, double size_stop, int n_iter)
{
    char str[SSM_STR_BUFFSIZE];

    double fitness = 0.0;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *simp = NULL;
    gsl_multimin_function minex_func;

    int iter = 0;
    int status;
    double size;

    gsl_vector *x = gsl_vector_alloc(nav->theta_all->length);
    gsl_vector *jump_sizes = gsl_vector_alloc(nav->theta_all->length);

    int i;
    for (i=0; i<nav->theta_all->length; i++) {
	gsl_vector_set(x, i, gsl_vector_get(theta, i));
	gsl_vector_set(jump_sizes, i, sqrt(gsl_matrix_get(var, i, i))); //note the sqrt !!
    }

    minex_func.n = nav->theta_all->length;
    minex_func.f = f_simplex;
    minex_func.params = params;

    simp = gsl_multimin_fminimizer_alloc(T, nav->theta_all->length);

    gsl_multimin_fminimizer_set(simp, &minex_func, x, jump_sizes);


    do {
	iter++;
	status = gsl_multimin_fminimizer_iterate(simp);
	if (status) break;
	size = gsl_multimin_fminimizer_size(simp);
	status = gsl_multimin_test_size(size, size_stop);

	fitness = - gsl_multimin_fminimizer_minimum(simp);

	if (nav->print & SSM_PRINT_LOG) {
	    if (status == GSL_SUCCESS) {
		ssm_print_log ("converged to maximum !");
	    }
	    sprintf(str, "%d\t logLike.: %12.5f\t size: %.14f", iter, fitness, size);
	    ssm_print_log(str);
	}

	if(nav->print & SSM_PRINT_TRACE){
	    ssm_print_trace(nav->trace, gsl_multimin_fminimizer_x(simp), nav, fitness, iter-1);
	}

    } while (status == GSL_CONTINUE && iter < n_iter);


    gsl_vector_memcpy(theta, gsl_multimin_fminimizer_x(simp));

    gsl_multimin_fminimizer_free(simp);
    gsl_vector_free(x);
    gsl_vector_free(jump_sizes);

    return fitness;
}
