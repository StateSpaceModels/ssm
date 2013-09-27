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

double ssm_log_likelihood(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav)
{
    int i;
    double loglike = 0.0;

    for(i=0; i< row->ts_nonan_length; i++){
	loglike += log(row->f_fitness[i](row->values[row->ts_nonan[i]], X, par, calc, t));
    }

    return loglike;
}


double ssm_sum_square(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav)
{
    int i;
    double ss = 0.0;

    for(i=0; i< row->ts_nonan_length; i++){
	ss += pow( row->values[row->ts_nonan[i]] - row->obs_mean[i](X, par, calc, t), 2 );
    }

    return ss;
}
