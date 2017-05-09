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
 ssm_err_code_t ssm_log_prob_proposal(double *log_proposal, ssm_theta_t *proposed, ssm_theta_t *theta, ssm_var_t *var, double sd_fac, ssm_nav_t *nav, int is_mvn)
 {

    int i;
    ssm_parameter_t *p;
    double p_tmp =0.0;
    double lp = 0.0;

    if (is_mvn) {
        p_tmp = ssm_dmvnorm(proposed->size, proposed, theta, var, sd_fac);
    }

    for(i=0; i<nav->theta_all->length; i++) {

        p = nav->theta_all->p[i];

        if (!is_mvn) {
            p_tmp = gsl_ran_gaussian_pdf((gsl_vector_get(proposed, p->offset_theta)-gsl_vector_get(theta, p->offset_theta)), sd_fac*sqrt(gsl_matrix_get(var, p->offset_theta, p->offset_theta)));
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

          p_tmp /= p->f_der_inv(gsl_vector_get(proposed, p->offset_theta));

        //check for numerical issues
          if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
            return SSM_ERR_PROPOSAL;
        }

        if (!is_mvn) {
            lp += log(p_tmp);
        }
    }

    if (is_mvn) {
        lp = log(p_tmp);
    }

    //check AGAIN for numerical issues (taking log could have created issues)
    if((isinf(lp)==1) || (isnan(lp)==1)) {
        return SSM_ERR_PROPOSAL;
    }

    *log_proposal = lp;

    return SSM_SUCCESS;
}


/**
 * log_prob_prior always compute a logged value (in the same way as
 * sanitize likelihood would have done) even if it fails (by that I
 * means it doesn't immediatly return on failure). This is usefull for
 * the --prior option.
 */
 ssm_err_code_t ssm_log_prob_prior(double *log_prior, ssm_theta_t *theta, ssm_nav_t *nav, ssm_fitness_t *fitness)
 {
    int i;
    ssm_parameter_t *p;
    int is_err = 0;
    double p_tmp = 0.0;
    double lp = 0.0;
    double back_transformed; //prior are on the natural scale, so we transform the parameter

    for(i=0; i<nav->theta_all->length; i++) {
        p = nav->theta_all->p[i];

        back_transformed = p->f_inv(gsl_vector_get(theta, p->offset_theta));
        p_tmp = p->f_prior(back_transformed);

        //check for numerical issues
        if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
            is_err =1;
            p_tmp = fitness->like_min;
        } else {
            p_tmp = (p_tmp <= fitness->like_min) ? fitness->like_min : p_tmp;
        }


        lp += log(p_tmp);
    }

    *log_prior = lp;

    return (is_err) ? SSM_ERR_PRIOR: SSM_SUCCESS;
}


/**
 * return accepted (SSM_SUCCESS) or rejected (SSM_MH_REJECTED) or
 * combination of prob errors and assign fitness->log_prior
 */
 ssm_err_code_t ssm_metropolis_hastings(ssm_fitness_t *fitness, double *alpha, ssm_theta_t *proposed, ssm_theta_t *theta, gsl_matrix *var, double sd_fac , ssm_nav_t *nav, ssm_calc_t *calc, int is_mvn)
 {
    double ran;
    ssm_err_code_t success = SSM_SUCCESS;
    double lproposal, lproposal_prev, lprior_prev;
    success |= ssm_log_prob_proposal(&lproposal,          proposed, theta,    var, sd_fac, nav,  is_mvn); /* q{ theta* | theta(i-1) }*/
    success |= ssm_log_prob_proposal(&lproposal_prev,     theta,    proposed, var, sd_fac, nav, is_mvn);  /* q{ theta(i-1) | theta* }*/
    success |= ssm_log_prob_prior   (&fitness->log_prior, proposed,                        nav, fitness); /* p{theta*} */
    success |= ssm_log_prob_prior   (&lprior_prev,        theta,                           nav, fitness); /* p{theta(i-1)} */


    if(success == SSM_SUCCESS) {

        // ( p{theta*}(y)  p{theta*} ) / ( p{theta(i-1)}(y) p{theta(i-1)} )  *  q{ theta(i-1) | theta* } / q{ theta* | theta(i-1) }
        *alpha = exp( (fitness->log_like - fitness->log_like_prev + lproposal_prev - lproposal + fitness->log_prior - lprior_prev) );
        ran = gsl_ran_flat(calc->randgsl, 0.0, 1.0);

        if(ran < *alpha) {
            return success; //accepted
        }
    } else {
        *alpha = 0.0;
        return success|SSM_MH_REJECT; //rejected
    }

    return SSM_MH_REJECT;
}

/**
 * return the empirical covariance matrix or the initial one
 * (depending on the iteration value and options) and the evaluated
 * tuning factor sd_fac
 */
 ssm_var_t *ssm_adapt_eps_var_sd_fac(double *sd_fac, ssm_adapt_t *a, ssm_var_t *var, ssm_nav_t *nav, int m)
 {
    // evaluate epsilon(m) = epsilon(m-1) * exp(a^(m-1) * (acceptance_rate(m-1) - 0.234))

    if ( (m > a->eps_switch) && ( m * a->ar < a->m_switch) ) {
        double ar = (a->flag_smooth) ? a->ar_smoothed : a->ar;
        a->eps *=  exp(pow(a->eps_a, (double) (m-1)) * (ar - 0.234));
    } else {  // after switching epsilon is set back to 1
        a->eps = 1.0;
    }

    a->eps = GSL_MIN(a->eps, a->eps_max);

    // evaluate tuning factor sd_fac = epsilon * 2.38/sqrt(n_to_be_estimated)
    *sd_fac = a->eps * 2.38/sqrt(nav->theta_all->length);

    return ((m * a->ar) >= a->m_switch) ? a->var_sampling: var;
}



/**
 * Compute acceptance rate using an average filter and smoothed
 * acceptance rate using a low pass filter (a.k.a exponential
 * smoothing or exponential moving average)
 */
void ssm_adapt_ar(ssm_adapt_t *a, int is_accepted, int m)
{
    a->ar_smoothed = (1.0 - a->alpha) * a->ar_smoothed + a->alpha * is_accepted;
    a->ar += ( ((double) is_accepted - a->ar) / ((double) (m + 1)) );
}



void ssm_theta_ran(ssm_theta_t *proposed, ssm_theta_t *theta, ssm_var_t *var, double sd_fac, ssm_calc_t *calc, ssm_nav_t *nav, int is_mvn)
{
    int i;

    if(is_mvn){
        ssm_rmvnorm(calc->randgsl, nav->theta_all->length, theta, var, sd_fac, proposed);
    } else { //iid gaussians
        for (i=0; i< nav->theta_all->length ; i++) {
            gsl_vector_set(proposed, i, gsl_vector_get(theta, i) + gsl_ran_gaussian(calc->randgsl, sd_fac*sqrt(gsl_matrix_get(var, i, i))));
        }
    }

}


int ssm_theta_copy(ssm_theta_t *dest, ssm_theta_t *src)
{
    return gsl_vector_memcpy(dest, src);
}

int ssm_par_copy(ssm_par_t *dest, ssm_par_t *src)
{
    return gsl_vector_memcpy(dest, src);
}




/**
 * The key is to understand that: X_resampled[j] = X[select[j]] so select
 * give the index of the resample ancestor...
 * The ancestor of particle j is select[j]
 *
 * With n index: X[n+1][j] = X[n][select[n][j]]
 *
 * Other caveat: D_J_p_X are in [N_DATA+1] ([0] contains the initial conditions)
 * select is in [N_DATA], times is in [N_DATA]
 */
 void ssm_sample_traj(ssm_X_t **D_X, ssm_X_t ***D_J_X, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness)
 {
    int j_sel; //the selected particle
    int n, nn, indn;

    double ran, cum_weights;

    ran=gsl_ran_flat(calc->randgsl, 0.0, 1.0);

    j_sel=0;
    cum_weights=fitness->weights[0];

    while (cum_weights < ran) {
        cum_weights += fitness->weights[++j_sel];
    }

    //print traj of ancestors of particle j_sel;

    //!!! below we assume that the last data point contain information'
    ssm_X_copy(D_X[data->n_obs], D_J_X[data->n_obs][j_sel]);

    // take ancestor of last data point
    j_sel = fitness->select[data->n_obs - 1][j_sel];

    //printing all ancesters up to previous observation time
    for(nn = (data->ind_nonan[data->n_obs_nonan-1]-1); nn > data->ind_nonan[data->n_obs_nonan-2]; nn--) {
        ssm_X_copy(D_X[nn + 1], D_J_X[nn + 1][j_sel]);
    }

    for(n = (data->n_obs_nonan-2); n >= 1; n--) {
        //indentifying index of the path that led to sampled particule
        indn = data->ind_nonan[n];
        ssm_X_copy(D_X[indn + 1], D_J_X[indn + 1][j_sel]);

        j_sel = fitness->select[indn][j_sel];

        //printing all ancesters up to previous observation time
        for(nn= (indn-1); nn > data->ind_nonan[n-1]; nn--) {
            ssm_X_copy(D_X[nn + 1], D_J_X[nn + 1][j_sel]);
        }
    }


    indn = data->ind_nonan[0];
    ssm_X_copy(D_X[nn + 1], D_J_X[ nn + 1 ][j_sel]);
    
    j_sel = fitness->select[indn][j_sel];

    for(nn=indn; nn>=-1; nn--) {
        ssm_X_copy(D_X[nn + 1], D_J_X[ nn + 1 ][j_sel]);
    }

}


// this is a function used for debugging a bug in retrieving particle genealogy
// it's similar to ssm_sample_traj but take the j_select as input
// j_select is returned by ssm_sample_traj_print2 and allow to compare the print and non-print version of this function
// we keep it in case we need to debugg it in the future
void ssm_sample_traj2(ssm_X_t **D_X, ssm_X_t ***D_J_X, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness,const int j_select)
{
    int j_sel; //the selected particle
    int n, nn, indn;

    double ran, cum_weights;

    ran=gsl_ran_flat(calc->randgsl, 0.0, 1.0);

    j_sel=0;
    cum_weights=fitness->weights[0];

    while (cum_weights < ran) {
        cum_weights += fitness->weights[++j_sel];
    }

    j_sel = j_select;

    // /*test*/
    // FILE *test_file = fopen("/Users/Tonton/work/projects/hev-modelling/ssm/SIR_ssm/pmcmc/X_sampled_0.csv", "a");
    // if (test_file == NULL)
    // {
    //     printf("Error opening file!\n");
    //     exit(1);
    // }

    // fprintf(test_file, "%i\n", j_sel);
    // fclose(test_file);
    // /*test*/

    //print traj of ancestors of particle j_sel;

    //!!! we assume that the last data point contain information'
    fprintf(stderr, "data->n_obs = %i\n", data->n_obs);
    ssm_X_copy(D_X[data->n_obs], D_J_X[data->n_obs][j_sel]);

    //printing all ancesters up to previous observation time
    fprintf(stderr, "data->n_obs_nonan = %i\n", data->n_obs_nonan);
    fprintf(stderr, "(data->ind_nonan[data->n_obs_nonan - 1] - 1) = %i\n", (data->ind_nonan[data->n_obs_nonan - 1] - 1));
    fprintf(stderr, "data->ind_nonan[data->n_obs_nonan-2] = %i\n", data->ind_nonan[data->n_obs_nonan-2]);

    // for(nn = 0; nn < data->n_obs_nonan; nn++){
    //     fprintf(stderr, "data->ind_nonan[%i] = %i\n", nn, data->ind_nonan[nn]); 
    // }

    j_sel = fitness->select[data->n_obs - 1][j_sel];

    for(nn = (data->ind_nonan[data->n_obs_nonan - 1] - 1); nn > data->ind_nonan[data->n_obs_nonan-2]; nn--) {
        fprintf(stderr, "nn = %i\n", nn);
        ssm_X_copy(D_X[nn+1], D_J_X[nn+1][j_sel]);
    }

    for(n = (data->n_obs_nonan-2); n >= 1; n--) {
        //indentifying index of the path that led to sampled particule
        indn = data->ind_nonan[n];
        ssm_X_copy(D_X[indn + 1], D_J_X[indn + 1][j_sel]);

        j_sel = fitness->select[indn][j_sel];

        //printing all ancesters up to previous observation time
        for(nn= (indn-1); nn > data->ind_nonan[n-1]; nn--) {
            ssm_X_copy(D_X[nn + 1], D_J_X[nn + 1][j_sel]);
        }
    }

    indn = data->ind_nonan[0];
    ssm_X_copy(D_X[nn + 1], D_J_X[ nn + 1 ][j_sel]);
    
    j_sel = fitness->select[indn][j_sel];

    for(nn=indn; nn>=-1; nn--) {
        ssm_X_copy(D_X[nn + 1], D_J_X[ nn + 1 ][j_sel]);
    }

    // j_sel = fitness->select[indn][j_sel];


    //TODO nn=-1 (for initial conditions)

}

