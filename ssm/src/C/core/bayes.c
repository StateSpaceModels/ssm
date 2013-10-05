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
 * compute the log of the prob of the proposal
 */
ssm_err_code_t ssm_log_prob_proposal(double *log_like, ssm_theta_t *proposed, ssm_theta_t *mean, ssm_var_t *var, double sd_fac, ssm_nav_t *nav, int is_mvn)
{

    int i, offset;
    ssm_parameter_t *p;
    double p_tmp =0.0;
    double Lp = 0.0;
    
    if (is_mvn) {
        p_tmp = ssm_dmvnorm(proposed->size, proposed, mean, var, sd_fac);
    }

    for(i=0; i<nav->theta_all->length; i++) {

	p = nav->theta_all->p[i];

	if (!is_mvn) {
	    p_tmp = gsl_ran_gaussian_pdf((gsl_vector_get(proposed, p->offset)-gsl_vector_get(mean, p->offset)), sd_fac*sqrt(gsl_matrix_get(var, p->offset, p->offset)));
	}

	/*
	  Change of variable formula:
	  Y = r(X) (r is assumed to be increasing) with inverse r_inv
	  X has f for density
	  Y has g for density

	  g(y) = f (r_inv (y)) * d r_inv(y)/dy
	  In our case, the proposal takes place on a transformed scale router.f. We know g(y) but we want f(r_inv(y))

	  we therefore divide g(y) by d r_inv(y)/dy.

	  Note that in the multivariate case, we need to
	  divide by the determinant of the Jacobian
	  matrix. However in our case the Jacobian is
	  diagonal so the determinant is the product of the
	  diagonal terms so everything generalizes nicely
	*/

	p_tmp /= p->f_inv_derivative(gsl_vector_get(proposed, p->offset));

	//check for numerical issues
	if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
	    return SSM_ERR_LIKE;
	}

	if (!is_mvn) {
	    Lp += log(p_tmp);
	}
    }

    if (is_mvn) {
	Lp = log(p_tmp);
    }

    //check AGAIN for numerical issues (taking log could have created issues)
    if((isinf(Lp)==1) || (isnan(Lp)==1)) {
	return SSM_ERR_LIKE;
    }
    
    *log_like = Lp;

    return SSM_SUCCESS;
}


/**
 * log_prob_prior always compute a logged value (in the same way as
 * sanitize likelihood would have done) even if it fails (by that I
 * means it doesn't immediatly return on failure). This is usefull for
 * the --prior option.
 */
ssm_err_code_t log_prob_prior(double *log_like, ssm_theta_t *mean, ssm_var_t *var, ssn_nav_t *nav, ssm_fitness_t *fitness)
{
    int i;
    ssm_parameter_t *p;
    int is_err = 0;
    double p_tmp = 0.0;
    double Lp = 0.0;
    double back_transformed; //prior are on the natural scale, so we transform the parameter

    for(i=0; i<nav->theta_all->length; i++) {

	p = nav->theta_all->p[i];

	back_transformed = p->f_inv(gsl_vector_get(mean, p->offset));
	p_tmp = p->prior(back_transformed);

	//check for numerical issues
	if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
	    is_err =1;    
	    p_tmp = fitness->like_min;
	} else {
	    p_tmp = (p_tmp <= fitness->like_min) ? fitness->like_min : p_tmp;
	}

	Lp += log(p_tmp); 
    }

    *log_like = Lp;

    return (is_err) ? SSM_ERR_LIKE: SSM_SUCCESS;
}


/**
 * return accepted (1) or rejected (0)
 */
int ssm_metropolis_hastings(double *alpha, ssm_theta_t *proposed, ssm_theta_t *mean, gsl_matrix *var, double sd_fac, ssm_fitness_t *fitness , ssm_nav_t *nav, ssm_calc_t *calc, int is_mvn)
{
    double ran;

    double lproposal_new, lproposal_prev, lprior_new, lprior_prev;
    ssm_err_code rc_proposal_new =  ssm_log_prob_proposal(&lproposal_new,  proposed, mean,     var, sd_fac, nav,  is_mvn); /* q{ theta* | theta(i-1) }*/
    ssm_err_code rc_proposal_prev = ssm_log_prob_proposal(&lproposal_prev, mean,     proposed, var, sd_fac, nav, is_mvn);  /* q{ theta(i-1) | theta* }*/
    ssm_err_code rc_prior_new  =    ssm_log_prob_prior(&lprior_new,        proposed,           var, nav, fitness);         /* p{theta*} */
    ssm_err_code rc_prior_prev =    ssm_log_prob_prior(&lprior_prev,       mean,               var, nav, fitness);         /* p{theta(i-1)} */

    if( (rc_proposal_new == SSM_SUCCESS) && (rc_proposal_prev == SSM_SUCCESS) && (rc_prior_new == SSM_SUCCESS) && (rc_prior_prev == SSM_SUCCESS) ) {

        // ( p{theta*}(y)  p{theta*} ) / ( p{theta(i-1)}(y) p{theta(i-1)} )  *  q{ theta(i-1) | theta* } / q{ theta* | theta(i-1) }
        *alpha = exp( (fitness->log_like_new - fitness->log_like_prev + lproposal_prev - lproposal_new + lprior_new - lprior_prev) );

        ran = gsl_ran_flat(p_calc->randgsl, 0.0, 1.0);

        if(ran < *alpha) {
            return 1; //accepted
        }
    } else {
        *alpha = 0.0;
        return 0; //rejected
    }

    return 0;
}
