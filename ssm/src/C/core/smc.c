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
 * Note that p_like->weights already contains the likelihood
 * @return the sucess status (sucess if some particles have a likelihood > LIKE_MIN)
 */
int ssm_weight(ssm_fitness_t *like, int n)
{

#if FLAG_WARNING
    char str[100];
#endif

    int j;

    double like_tot_n = 0.0;
    int nfailure_n = 0;
    int success = 1;

    like->ess_n = 0.0;

    for(j=0; j < like->J ; j++) {
        /*compute first part of weights (non divided by sum likelihood)*/
        if (like->weights[j] <= pow(LIKE_MIN, N_TS)) {
            like->weights[j] = 0.0;
            nfailure_n += 1;
        } else {
            like_tot_n += like->weights[j]; //note that like_tot_n contains only like of part having a like>LIKE_MIN
            like->ess_n += like->weights[j]*like->weights[j]; //first part of ess computation (sum of square)
        }
    } /*end of for on j*/

    /*compute second part of weights (divided by sum likelihood)*/
    if(nfailure_n == like->J) {
        success = 0;
        like->n_all_fail += 1;
#if FLAG_WARNING
        sprintf(str,"warning: nfailure = %d, at n=%d we keep all particles and assign equal weights", nfailure_n , n);
        print_warning(str);
#endif
        like->log_like_n = LOG_LIKE_MIN*N_TS;

        double invJ=1.0/ ((double) like->J);
        for(j=0 ; j < like->J ; j++) {
            like->weights[j]= invJ;
            like->select[n][j]= j;
        }
        like->ess_n = 0.0;

    } else {
        for(j=0 ; j < like->J ; j++) {
            like->weights[j] /= like_tot_n;
        }

        like->log_like_n = log(like_tot_n / ((double) like->J));
        like->ess_n = (like_tot_n*like_tot_n)/like->ess_n;
    }

    like->log_like += like->log_like_n;

    return success;
}


/**
 *   Systematic sampling.  Systematic sampling is faster than
 *   multinomial sampling and introduces less monte carlo variability
 */
void ssm_systematic_sampling(ssm_fitness_t *like, ssm_calc_t *calc, int n)
{
    unsigned int *select = like->select[n];
    double *prob = like->weights;

    int i,j;

    double ran;
    double inc = 1.0/((double) like->J);

    ran = gsl_ran_flat(calc->randgsl, 0.0, inc);
    i = 0;
    double weight_cum = prob[0];

    for(j=0; j < like->J; j++) {
        while(ran > weight_cum) {
            i++;
            weight_cum += prob[i];
        }
        select[j] = i;
        ran += inc;
    }
}


void ssm_resample_X(ssm_fitness_t *like, struct s_X ***J_p_X, struct s_X ***J_p_X_tmp, int n)
{
    int k, j;

    unsigned int *select = like->select[n];

    for(j=0;j<like->J;j++) {

        (*J_p_X_tmp)[j]->dt = (*J_p_X)[select[j]]->dt;

	gsl_vector_memcpy((*J_p_X_tmp)[j]->proj, (*J_p_X)[select[j]]->proj);	
    }

    swap_X(J_p_X, J_p_X_tmp);
}


void swap_X(struct s_X ***X, struct s_X ***tmp_X)
{
    /* swap array of pointers of pointers to struct s_X */

    struct s_X **tmp;
    tmp=*tmp_X;
    *tmp_X=*X;
    *X=tmp;
}
