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

void ssm_X_copy(ssm_X_t *dest, ssm_X_t *src){
    gsl_vector_memcpy(dest->proj, src->proj);
    dest->dt = src->dt;
}

void ssm_X_reset_inc(ssm_X_t *X, ssm_data_row_t *row)
{
    int i;

    for(i=0; i<row->states_reset_length; i++){
	gsl_vector_set(X->proj, row->states_reset[i]->offset, 0.0);
    }
}

/**
 * Modified version of gsl_ran_multinomial to avoid a loop. We avoid
 * to recompute the total sum of p (called norm in GSL) as it will
 * always be 1.0 with ssm (no rounding error by construction)
 */
void ssm_ran_multinomial (const gsl_rng * r, const size_t K, unsigned int N, const double p[], unsigned int n[])
{
    size_t k;
    double sum_p = 0.0;

    unsigned int sum_n = 0;

    for (k = 0; k < K; k++) {
        if (p[k] > 0.0) {
            n[k] = gsl_ran_binomial (r, p[k] / (1.0 - sum_p), N - sum_n);
        }
        else {
            n[k] = 0;
        }

        sum_p += p[k];
        sum_n += n[k];
    }
}
