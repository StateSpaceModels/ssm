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
 * Computation of the EKF gain kt for observation data_t_ts and obs
 * jacobian ht, given estimate xk_t_ts and current covariance Ct
 */
ssm_err_code ssm_kalman_gain_computation(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav){
    
    int i,j;
    int cum_status = 0;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;

    // sub-matrices and sub-vectors of working variables are extracted as not all tseries are observed
    gsl_vector_view *pred_error = gsl_vector_subvector(calc->_pred_error,0,row->ts_nonan_length);
    gsl_matrix_view *St = gsl_matrix_submatrix(calc->_St,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view *Stm1 = gsl_matrix_submatrix(calc->_Stm1,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view *Rt = gsl_matrix_submatrix(calc->_Rt,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view *Tmp = gsl_matrix_submatrix(calc->_Tmp_N_SV_N_TS,0,0,m,row->ts_nonan_length);
    gsl_matrix_view *Ht = gsl_matrix_subvector(calc->_Ht,0,0,m,row->ts_nonan_length);
    gsl_matrix_view *Kt = gsl_matrix_submatrix(calc->_Kt,0,0,m,row->ts_nonan_length);
    gsl_matrix_view *Ct =  gsl_matrix_const_view_array(&X[m], m, m);


    // fill Ht and Rt
    eval_Ht(X, row, par, nav, calc, t);
    for(i=0; i< row->ts_nonan_length; i++){
	for(j=0; j< row->ts_nonan_length; j++){
	    if (i==j){
		gsl_matrix_set(Rt,i,j) = row->observed[i]->obs_var(X, par, calc, t); 
	    } else {
		gsl_matrix_set(Rt,i,j) = 0;
	    }
	}
    }

    // pred_error = double data_t_ts - xk_t_ts
    for(i=0; i< row->ts_nonan_length; i++){
	gsl_vector_set(pred_error,i) = row->values[i] - row->observed[i]->obs_mean(X, par, calc, t);
    }


    // sc_st = Ht' * Ct * Ht + sc_rt
    /*
     * here ht is a column vector to fit gsl standards,
     * rather than a row vector as it should be in theory,
     * which explains the slightly different formula we use
     */

    // workn = Ct*Ht
    cum_status |= gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ct, Ht, 0.0, Tmp); 

    // sc_st = Ht' * workn;
    cum_status |= gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ht, Tmp, 0.0, St); 

    // sc_st = sc_st + sc_rt ;
    cum_status |=  gsl_matrix_add(St,Rt);


    // Kt = Ct * Ht * sc_st^-1
    cum_status |= gsl_linalg_LU_decomp(St, p, signum); // inversion requires LU decomposition
    cum_status |= gsl_linalg_LU_invert(St,p,Stm1);
    cum_status |= gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ht, Stm1, 0.0, Tmp);
    cum_status |= gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ct, Tmp, 0.0, Kt);

    return cum_status ? SSM_ERR_KAL : SSM_SUCCESS;

}




ssm_err_code ssm_kalman_update(ssm_X_t *X, ssm_row_t *row, double t, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *like){
    
    status = ssm_kalman_gain_computation(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav);    
    int cum_status = 0;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    gsl_vector_view *pred_error = gsl_matrix_subvector(calc->_pred_error,0,row->ts_nonan_length);
    gsl_matrix_view *Kt = gsl_matrix_submatrix(calc->_Kt,0,0,m,row->ts_nonan_length);
    gsl_matrix_view *Tmp = gsl_matrix_submatrix(calc->_Tmp_N_TS_N_SV,0,0,row->ts_nonan_length,m);
    gsl_vector_view *X_sv = gsl_matrix_subvector(X->proj,0,m);
    gsl_matrix_view *Ct =  gsl_matrix_const_view_array(&X[m], m, m);

    //////////////////
    // state update //
    //////////////////
    // X_sv += Kt * pred_error
    cum_status |= gsl_blas_dgemv(CblasNoTrans,1.0,Kt,pred_error,1.0,X_sv);

    ///////////////////////
    // covariance update //
    ///////////////////////
    // Ct = Ct - Kt * Ht' * Ct
    cum_status |= gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ht, Ct, 0.0, Tmp);
    cum_status |= gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, Kt, Tmp, 0.0, Ct);

    return status | (cum_status ? SSM_ERR_KAL : SSM_SUCCESS);
}
