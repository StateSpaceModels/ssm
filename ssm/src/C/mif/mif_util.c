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

double ssm_mif_cooling(ssm_options_t *opts, int m)
{
    return pow(opts->a, (double) (m-1));
}

/**
 * if x is a prameter and y an initial condition and n the number of data points:
 *  rescale cov(x,y) by 1/sqrt(n)
 *  rescale cov(x,x) by (1/sqrt(n))*(1/sqrt(n))
 *  do not rescale cov(y,y)
 */
void ssm_mif_scale_var(ssm_var_t *var, ssm_data_t *data, ssm_nav_t *nav)
{
    int i, j;
    int offset_i, offset_j;

    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    ssm_it_parameters_t *fls = nav->theta_icsv_icdiff;       //parameters fitted with fixed lag smoothing (fls)

    double inv = 1.0/((double) data->n_obs);
    double sqrt_inv = 1.0/sqrt((double) data->n_obs);

    //mif, mif terms: rescale by inv
    for(i=0; i<mif->length; i++){
        offset_i = mif->p[i]->offset_theta;
        for(j=0; j<mif->length; j++){
            offset_j = mif->p[j]->offset_theta;
            gsl_matrix_set(var, offset_i, offset_j, gsl_matrix_get(var, offset_i, offset_j) * inv);
        }
    }

    //mif, fls and fls, mif terms: rescale by sqrt_inv
    for(i=0; i<mif->length; i++){
        offset_i = mif->p[i]->offset_theta;
        for(j=0; j<fls->length; j++){
            offset_j = fls->p[j]->offset_theta;
            gsl_matrix_set(var, offset_i, offset_j, gsl_matrix_get(var, offset_i, offset_j) * sqrt_inv);
            gsl_matrix_set(var, offset_j, offset_i, gsl_matrix_get(var, offset_j, offset_i) * sqrt_inv);
        }
    }

}



/**
 * Multiply the loglikelihood of particle j by prod_i prior(theta_i)^(1/n_obs)
 */
void ssm_mif_patch_like_prior(double *like, ssm_fitness_t *fitness, ssm_theta_t **J_theta, ssm_data_t *data, ssm_nav_t *nav, const int n, const int lag)
{
    int i, j;

    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    ssm_it_parameters_t *fls = nav->theta_icsv_icdiff;       //parameters fitted with fixed lag smoothing (fls)

    ssm_parameter_t *p;
    double log_like_prior_i;
    double log_like;
    double inv_n_obs = 1.0/ ((double) data->n_obs);
    double inv_lag = 1.0/ ((double) lag);

    for(j=0; j<fitness->J; j++) {
        if(like[j]){
            log_like = log(like[j]);

            // likelihood is multiplied by prior(theta_j)^(1/n_obs) for parameters fitted with MIF (as opposed to fixed lag smoothing)
            for(i=0; i<mif->length; i++) {
                p = mif->p[i];
                log_like_prior_i = log(p->prior( p->f_inv(gsl_vector_get(J_theta[j], p->offset_theta)) ));
                log_like += inv_n_obs * log_like_prior_i;
            }

            // likelihood is multiplied by prior(theta_j)^(1/lag) for parameters fitted with fixed lag smoothing
            if(n<lag){
                for(i=0; i<fls->length; i++) {
                    p = fls->p[i];
                    log_like_prior_i = log(p->prior( p->f_inv(gsl_vector_get(J_theta[j], p->offset_theta)) ));
                    log_like += inv_lag * log_like_prior_i;
                }
            }

            //check for numerical issues and convert back to likelihood
            if( (isnan(log_like) == 1) || (isinf(log_like) == 1) || (log_like > 1.0) ) {
                like[j] = 0.0;
            } else {
                like[j] = exp(log_like);
            }
        }
    }
}


/**
 * Compute filtered mean and prediction var of particles at time
 * n. We take weighted averages with "weights" for the filtered mean
 * (in order to reduce monte-carlo variability) and use a numericaly
 * stable online algo for the variance.
 */
void ssm_mif_mean_var_theta_theoretical(double *theta_bart, double *theta_Vt, ssm_theta_t **J_theta, ssm_var_t *var, ssm_fitness_t *fitness, ssm_nav_t *nav, double var_fac)
{
    int i, j;
    double kn, M2, avg, delta; //for variance computations
    int offset;

    ssm_it_parameters_t *it;
    if (nav->print & SSM_PRINT_DIAG){
        it= nav->theta_all;
    } else {
        it= nav->theta_no_icsv_no_icdiff; //only this one is truely needed
    }

    for(i=0; i<it->length; i++) {
        offset = it->p[i]->offset_theta;
        theta_bart[offset]=0.0;

        kn=0.0;
        avg=0.0;
        M2=0.0;

        for(j=0 ; j<fitness->J ; j++) {
            //variance computation
            kn += 1.0;
            delta = gsl_vector_get(J_theta[j], offset) - avg;
            avg += (delta / kn);
            M2 += delta*(gsl_vector_get(J_theta[j], offset) - avg);

            //weighted average for filtered mean
            theta_bart[offset] += fitness->weights[j]*gsl_vector_get(J_theta[j], offset);
        }
        theta_Vt[offset] = M2/(kn -1.0);

        if( (theta_Vt[offset]<0.0) || (isinf(theta_Vt[offset])==1) || (isnan(theta_Vt[offset])==1)) {
            theta_Vt[offset]=0.0;

            if (!(nav->print & SSM_QUIET)) {
                ssm_print_warning("error in prediction variance computation");
            }
        }

        /*we add theoretical variance corresponding to
          mutation of theta to reduce Monte Carlo variability
          AND ensure that Vt is > 0.0 (so that the crappy
          Ionides formulae don't crash (even if it will even
          with that)')*/
        //TODO check that this is valid in MVN case
        theta_Vt[offset] += var_fac*gsl_matrix_get(var, offset, offset);
    }
}


void ssm_mif_resample_and_mutate_theta(ssm_fitness_t *fitness, ssm_theta_t **J_theta, ssm_theta_t **J_theta_tmp, ssm_var_t *var, ssm_calc_t **calc, ssm_nav_t *nav, double sd_fac, int n)
{
    int i, j, offset;

    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    ssm_it_parameters_t *fls = nav->theta_icsv_icdiff;       //parameters fitted with fixed lag smoothing (fls)
    unsigned int *select = fitness->select[n];

    //resample
    for(j=0; j<fitness->J; j++) {
        ssm_theta_copy(J_theta_tmp[j], J_theta[select[j]]);
    }

    for(j=0; j<fitness->J; j++) {
        //resample and mutate
        for(i=0; i<mif->length; i++) {
            offset = mif->p[i]->offset_theta;
            gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset) + gsl_ran_gaussian(calc[0]->randgsl, sd_fac*sqrt(gsl_matrix_get(var, offset, offset))));
        }

        //resample only
        for(i=0; i<fls->length; i++) {
            offset = fls->p[i]->offset_theta;
            gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset));
        }
    }
}



void ssm_mif_fixed_lag_smoothing(ssm_theta_t *mle, ssm_theta_t **J_theta, ssm_fitness_t *fitness, ssm_nav_t *nav)
{
    int i, j;
    ssm_it_parameters_t *fls = nav->theta_icsv_icdiff;

    for(i=0; i<fls->length; i++) {
        int offset = fls->p[i]->offset_theta;
        gsl_vector_set(mle, offset, 0.0);
        for(j=0; j<fitness->J; j++) {
            gsl_vector_set(mle, offset, gsl_vector_get(mle, offset) + gsl_vector_get(J_theta[j], offset) * fitness->weights[j]);
        }
    }
}


/**
 * Update mle in a numericaly stable way for the first
 * iterations (as suggested in Ionides et al. 2006)
 *
 * NOTE: D_theta_bart is in [data_length+1] so everything has to be
 * shiftted by 1
 */
void ssm_mif_update_average(ssm_theta_t *mle, double **D_theta_bart, ssm_data_t *data, ssm_nav_t *nav)
{
    int i, n, offset;
    double tmp;
    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff;

    for(i=0; i<mif->length; i++) {
        offset = mif->p[i]->offset_theta;
        tmp = 0.0;
        for(n=0; n<data->n_obs_nonan; n++){
            tmp += D_theta_bart[data->ind_nonan[n] + 1][offset];
        }
        gsl_vector_set(mle, offset, tmp / ((double) data->n_obs_nonan) );
    }
}


/**
 * The MIF update formulae Ionides et al 2006 PNAS (doesn't work
 * although the authors said it should)
 *
 * theta_{m+1} = theta_m + V(t1) * sum_n{1...n} [ 1/V(t_n) (theta(t_n) - theta(t_{n-1})) ]
 *
 * NOTE: D_theta_bart is in [data_length+1] so everything has to be
 * shiftted by 1
 */
void ssm_mif_update_ionides(ssm_theta_t *mle, ssm_var_t *var, double **D_theta_bart, double **D_theta_Vt, ssm_data_t *data, ssm_nav_t *nav, ssm_options_t *opts, double cooling)
{
    int i,  n, nn, nnp1, offset;
    double tmp;
    ssm_it_parameters_t *mif = nav->theta_no_icsv_no_icdiff;

    for(i=0; i<mif->length; i++) {
        offset = mif->p[i]->offset_theta;

        //from initial condition (before first data point) to first data point
        nnp1 = data->ind_nonan[0]; //first entry with data
        tmp = ( (D_theta_bart[nnp1 + 1][offset] - D_theta_bart[0][offset]) / D_theta_Vt[nnp1+1][offset] );

        //from data point to data point
        for(n=1; n< data->n_obs_nonan; n++){
            nn = data->ind_nonan[n-1];  //previous entry with data
            nnp1 = data->ind_nonan[n];  //current entry with data
            tmp += (D_theta_bart[nnp1+1][offset] - D_theta_bart[nn+1][offset]) / D_theta_Vt[nnp1+1][offset];
        }

        gsl_vector_set(mle,
                       offset,
                       (gsl_vector_get(mle, offset) + ((data->rows[0]->time + opts->b*opts->b)*cooling*cooling*gsl_matrix_get(var, offset, offset)*tmp) ) );  //Cf Ionides et al PNAS 2006
    }
}


void ssm_mif_print_mean_var_theoretical_ess(FILE *stream, double *theta_bart, double *theta_Vt, ssm_fitness_t *fitness, ssm_nav_t *nav , ssm_row_t *row, int m)
{

    int i, offset;
    char key[SSM_STR_BUFFSIZE];

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(m));
    json_object_set_new(jout, "date", json_string(row->date));

    for(i=0; i<nav->theta_all->length; i++){
        offset = nav->theta_all->p[i]->offset_theta;

        json_object_set_new(jout, nav->theta_all->p[i]->name, json_real(theta_bart[offset]));
        snprintf(key, SSM_STR_BUFFSIZE, "var_%s", nav->theta_all->p[i]->name);
        json_object_set_new(jout, key, json_real(theta_Vt[offset]));

        json_object_set_new(jout, "ess", isnan(fitness->ess_n) ? json_null() : json_real(fitness->ess_n));
    }

    ssm_json_dumpf(stream, "mif", jout);
}
