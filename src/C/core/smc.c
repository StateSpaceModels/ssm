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
 * Computes the weight of the particles.
 * Note that p_fitness->weights already contains the likelihood
 * @return the sucess status (sucess if some particles have a likelihood > LIKE_MIN)
 */
int ssm_weight(ssm_fitness_t *fitness, ssm_row_t *row, ssm_nav_t *nav, int n)
{
    char str[SSM_STR_BUFFSIZE];

    int j;

    double like_tot_n = 0.0;
    int nfailure_n = 0;
    int success = 1;

    fitness->ess_n = 0.0;

    for(j=0; j < fitness->J ; j++) {
        /*compute first part of weights (non divided by sum likelihood)*/      
        if (fitness->weights[j] <= pow(fitness->like_min, row->ts_nonan_length)) {
            fitness->weights[j] = 0.0;
            nfailure_n += 1;
        } else {
            like_tot_n += fitness->weights[j]; //note that like_tot_n contains only like of part having a like>like_min
            fitness->ess_n += fitness->weights[j]*fitness->weights[j]; //first part of ess computation (sum of square)
        }
    } /*end of for on j*/

    /*compute second part of weights (divided by sum likelihood)*/
    if(nfailure_n == fitness->J) {
        success = 0;
        fitness->n_all_fail += 1;
        if (nav->print & SSM_PRINT_WARNING) {
            sprintf(str,"nfailure = %d, at n=%d we keep all particles and assign equal weights", nfailure_n, n);
            ssm_print_warning(str);
        }
        fitness->log_like_n = fitness->log_like_min * row->ts_nonan_length;

        double invJ=1.0/ ((double) fitness->J);
        for(j=0 ; j < fitness->J ; j++) {
            fitness->weights[j]= invJ;
            fitness->select[n][j]= j;
        }
        fitness->ess_n = 0.0;

    } else {
        for(j=0 ; j < fitness->J ; j++) {
            fitness->weights[j] /= like_tot_n;
        }

        fitness->log_like_n = log(like_tot_n / ((double) fitness->J));
        fitness->ess_n = (like_tot_n*like_tot_n)/fitness->ess_n;
    }

    fitness->log_like += fitness->log_like_n;

    return success;
}


/**
 *   Systematic sampling.  Systematic sampling is faster than
 *   multinomial sampling and introduces less monte carlo variability
 */
void ssm_systematic_sampling(ssm_fitness_t *fitness, ssm_calc_t *calc, int n)
{
    unsigned int *select = fitness->select[n];
    double *prob = fitness->weights;

    int i,j;

    double ran;
    double inc = 1.0/((double) fitness->J);

    ran = gsl_ran_flat(calc->randgsl, 0.0, inc);
    i = 0;
    double weight_cum = prob[0];

    for(j=0; j < fitness->J; j++) {
        while(ran > weight_cum) {
            i++;
            weight_cum += prob[i];
        }
        select[j] = i;
        ran += inc;
    }
}


void ssm_resample_X(ssm_fitness_t *fitness, ssm_X_t ***J_p_X, ssm_X_t ***J_p_X_tmp, int n)
{
    int j;

    unsigned int *select = fitness->select[n];

    for(j=0;j<fitness->J;j++) {
        (*J_p_X_tmp)[j]->dt = (*J_p_X)[select[j]]->dt;
        memcpy((*J_p_X_tmp)[j]->proj, (*J_p_X)[select[j]]->proj, (*J_p_X)[select[j]]->length * sizeof(double));
    }

    ssm_swap_X(J_p_X, J_p_X_tmp);
}

/**
 * swap array of pointers of pointers to ssm_X_t
 */
void ssm_swap_X(ssm_X_t ***X, ssm_X_t ***tmp_X)
{
    ssm_X_t **tmp;
    tmp=*tmp_X;
    *tmp_X=*X;
    *X=tmp;
}
