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
 * Adapted from: Multivariate Normal density function and random
 * number generator Using GSL from Ralph dos Santos Silva
 * Copyright (C) 2006

 * multivariate normal distribution random number generator
 *
 * @param n      dimension of the random vetor
 * @param mean   vector of means of size n
 * @param var    variance matrix of dimension n x n
 * @param result output variable with a sigle random vector normal distribution generation
 */
int ssm_rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, double sd_fac, gsl_vector *result)
{

    int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(work,var);
    //scale var with sd_fac^2
    gsl_matrix_scale(work, sd_fac*sd_fac);

    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
        gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
    gsl_vector_add(result,mean);

    gsl_matrix_free(work);

    return 0;
}


/**
 * Adapted from: Multivariate Normal density function and random
 * number generator Using GSL from Ralph dos Santos Silva
 * Copyright (C) 2006
 *
 * multivariate normal density function
 *
 * @param n	dimension of the random vetor
 * @param mean	vector of means of size n
 * @param var	variance matrix of dimension n x n
 */
double ssm_dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, double sd_fac)
{
    int s;
    double ax,ay;
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n),
        *winv = gsl_matrix_alloc(n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    gsl_matrix_memcpy( work, var );
    //scale var with sd_fac^2
    gsl_matrix_scale(work, sd_fac*sd_fac);

    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );

    xm = gsl_vector_alloc(n);
    gsl_vector_memcpy( xm, x);
    gsl_vector_sub( xm, mean );
    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);
    ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

    return ay;
}



/**
 * evaluate the empirical covariance matrix using a stable one-pass
 * algorithm: see
 * http://en.wikipedia.org/wiki/Algorithms%5Ffor%5Fcalculating%5Fvariance#Covariance
 *
 * @param x       current best estimate of the parameters;
 * @param  m the current mcmc iteration index
 */
void ssm_adapt_var(ssm_adapt_t *adapt, ssm_theta_t *x, int m)
{
    int i, k;

    double *x_bar = adapt->mean_sampling;   
    gsl_matrix *cov = adapt->var_sampling;
    double dm = (double) m;
    double val;

    /* lower triangle and diagonal */
    for (i=0; i < x->size; i++) {
        for (k=0; k <= i; k++) {
	    //Cov_n = C_n/n
	    //C_n = sum_{i=1,n} (x_i-x_bar_n)(y_i-y_bar_n)
	    //C_n = C_{n-1} + (n-1)/n(x_n -x_bar_{n-1})(y_n -y_bar_{n-1})
	    //val = (n-1)/n(x_n -x_bar_{n-1})(y_n -y_bar_{n-1})
            val = ((dm - 1.0) / dm) * (gsl_vector_get(x, i) - x_bar[i]) * (gsl_vector_get(x, k) - x_bar[k]);

	    //C_n = C_{n-1} + val
	    //NOTE: we track the covariance and not C_n. However C_n =  covariance * n  so C_{n-1} = cov_{n-1} * (n-1) hence the fomula below
            gsl_matrix_set(cov, i, k, (((dm-1.0)*gsl_matrix_get(cov, i, k)) + val ) / dm);
        }
    }
    
    /* fill upper triangle */
    for (i=0; i < x->size; i++) {
        for (k=0; k < i; k++) {
            val = gsl_matrix_get(cov, i, k);
            gsl_matrix_set(cov, k, i, val);
        }
    }

    /* update x_bar */
    for (i=0; i < x->size; i++) {
        x_bar[i] += (gsl_vector_get(x, i) - x_bar[i]) / dm;
    }

}
