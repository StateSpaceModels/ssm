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
 *  checks for numerical issues.  Also, ensures a uniform likelihood
 *  scale by making sure that everything below fitness->log_like_min
 *  is fitness->log_like_min*row->ts_nonan_length.
 */
 double ssm_sanitize_log_likelihood(double log_like, ssm_row_t *row, ssm_fitness_t *fitness, ssm_nav_t *nav)
 {
    if ((isinf(log_like)==1) || (isnan(log_like)==1) ) { //error
        if (nav->print & SSM_PRINT_WARNING) {
            ssm_print_warning("error likelihood computation");
            // fprintf(stderr, "at time %u\n", row->time);
        }
        return fitness->log_like_min * row->ts_nonan_length;
    } else {
        return GSL_MAX(fitness->log_like_min * row->ts_nonan_length, log_like);
    }
}


double ssm_log_likelihood(ssm_row_t *row, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *fitness)
{
    int i;
    double loglike = 0.0;
    double t = row->time;

    for(i=0; i< row->ts_nonan_length; i++){
        loglike += log(row->observed[i]->f_likelihood(row->values[i], X, par, calc, t));
    }

    return ssm_sanitize_log_likelihood(loglike, row, fitness, nav);
}


double ssm_sum_square(ssm_row_t *row, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *fitness)
{
    int i;
    double ss = 0.0;
    double t = row->time;

    for(i=0; i< row->ts_nonan_length; i++){
        ss += pow( row->values[i] - row->observed[i]->f_obs_mean(X, par, calc, t), 2 );
    }

    return ss;
}



void ssm_aic(ssm_fitness_t *fitness, ssm_nav_t *nav, double log_like)
{
    int p = nav->theta_all->length;
    fitness->summary_log_likelihood = log_like;
    fitness->AIC = 2*p-2*log_like;
    fitness->AICc = fitness->AIC + 2*p*(p+1)/ ((double) (fitness->n - p -1));
    fitness->DIC = -2*log_like; //deviance in this case
}


void ssm_dic_init(ssm_fitness_t *fitness, double log_like, double log_prior)
{
    fitness->summary_log_likelihood = log_like;
    fitness->summary_log_ltp = log_like + log_prior;
    fitness->_min_deviance = -2*log_like;
    fitness->_deviance_cum = -2*log_like;
}


void ssm_dic_update(ssm_fitness_t *fitness, double log_like, double log_prior)
{
    if( log_like > fitness->summary_log_likelihood){
       fitness->summary_log_likelihood = log_like;
   }
   if( (log_like+log_prior) > fitness->summary_log_ltp){
       fitness->summary_log_ltp = (log_like+log_prior);
   }

   double deviance = -2*log_like;
   fitness->_deviance_cum += deviance;
   if( deviance < fitness->_min_deviance){
       fitness->_min_deviance = deviance;
   }
}

void ssm_dic_end(ssm_fitness_t *fitness, ssm_nav_t *nav, int m)
{
    //populate AIC fields (based on the best log_likelihood found in the trace)
    ssm_aic(fitness, nav, fitness->summary_log_likelihood);
    //overwrite DIC with correct value
    double mean_dev = fitness->_deviance_cum / ((double) m);
    fitness->DIC = 2*mean_dev - fitness->_min_deviance;
}
