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
 * stepping function for the Extended Kalman Filter
 */
void step_kalman(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    double *X = p_X->proj;
    double dt = p_X->dt;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    int i, j;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;

    gsl_matrix *Q = calc->_Q;
    gsl_matrix_view Ct   = gsl_matrix_view_array(&X[m], m, m);
    gsl_matrix *FtCt = calc->_FtCt;
    gsl_matrix *Ft = calc->_Ft;

    eval_Q(X, par, nav, calc, t);
    eval_jac(X, par, nav, calc, t);
    
    /*State variables propagation*/
    /*copy and paste ode*/

    /*Covariance propagation*/
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ft, &Ct.matrix, 0.0, FtCt);
    
    for(i=0; i< m; i++){
	for(c=0; j< m; j++){
	    gsl_matrix_set(&Ct, 
			   i,
			   j, 
			   gsl_matrix_get(FtCt, i, j) + gsl_matrix_get(FtCt, j, i) + gsl_matrix_get(Q, i, j));
	}
    }
}


