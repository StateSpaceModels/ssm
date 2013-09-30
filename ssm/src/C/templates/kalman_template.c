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
 * stepping functions for the Extended Kalman Filter
 */
void step_kalman(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    double *X = p_X->proj;
    double dt = p_X->dt;
    int n = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    int i, j;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;

    gsl_matrix *Q = calc->Q;
    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&X[n], n, n);
    gsl_matrix *FtCt = nav->FtCt;

    eval_Q(Q, X, par, nav, calc, t);
    eval_jac(X, par, nav, calc, t);
    
    /*State variables propagation*/
    /*copy and paste ode*/

    /*Covariance propagation*/
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, X->Ft, &Ct.matrix, 0.0, FtCt);
    for(i=0; i< n; i++){
	for(c=0; j< n; j++){
	    gsl_matrix_set(&Ct, 
			   i,
			   j, 
			   gsl_matrix_get(FtCt, i, j) + gsl_matrix_get(FtCt, j, i) + gsl_matrix_get(Q, i, j));
	}
    }
}


