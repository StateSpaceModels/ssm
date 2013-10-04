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
ssm_err_code_t ssm_kalman_gain_computation(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav){

    int i, j, status;
    ssm_err_code_t cum_status = SSM_SUCCESS;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;

    // sub-matrices and sub-vectors of working variables are extracted as not all tseries are observed
    gsl_vector_view pred_error = gsl_vector_subvector(&calc->_pred_error[0],0,row->ts_nonan_length);
    gsl_matrix_view St = gsl_matrix_submatrix(calc->_St,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view Stm1 = gsl_matrix_submatrix(calc->_Stm1,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view Rt = gsl_matrix_submatrix(calc->_Rt,0,0,row->ts_nonan_length,row->ts_nonan_length);
    gsl_matrix_view Tmp = gsl_matrix_submatrix(calc->_Tmp_N_SV_N_TS,0,0,m,row->ts_nonan_length);
    gsl_matrix_view Ht = gsl_matrix_submatrix(calc->_Ht,0,0,m,row->ts_nonan_length);
    gsl_matrix_view Kt = gsl_matrix_submatrix(calc->_Kt,0,0,m,row->ts_nonan_length);
    gsl_matrix_view Ct =  gsl_matrix_view_array(&X->proj[m], m, m);
    gsl_permutation *p  = gsl_permutation_alloc(m);

    // fill Ht and Rt
    ssm_eval_Ht(X, row, t, par, nav, calc);
    for(i=0; i< row->ts_nonan_length; i++){
        for(j=0; j< row->ts_nonan_length; j++){
            if (i==j){
                gsl_matrix_set(&Rt.matrix,i,j,row->observed[i]->f_obs_var(X, par, calc, t));
            } else {
                gsl_matrix_set(&Rt.matrix,i,j,0);
            }
        }
    }

    // pred_error = double data_t_ts - xk_t_ts
    for(i=0; i< row->ts_nonan_length; i++){
        gsl_vector_set(&pred_error.vector,i,row->values[i] - row->observed[i]->f_obs_mean(X, par, calc, t));
    }


    // sc_st = Ht' * Ct * Ht + sc_rt
    /*
     * here ht is a column vector to fit gsl standards,
     * rather than a row vector as it should be in theory,
     * which explains the slightly different formula we use
     */

    // workn = Ct*Ht
    status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Ct.matrix, &Ht.matrix, 0.0, &Tmp.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    // sc_st = Ht' * workn;
    status = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Ht.matrix, &Tmp.matrix, 0.0, &St.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    // sc_st = sc_st + sc_rt ;
    status  =  gsl_matrix_add(&St.matrix,&Rt.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    // Kt = Ct * Ht * sc_st^-1
    status = gsl_linalg_LU_decomp(&St.matrix, p, &i); // inversion requires LU decomposition
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    status = gsl_linalg_LU_invert(&St.matrix,p,&Stm1.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Ht.matrix, &Stm1.matrix, 0.0, &Tmp.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Ct.matrix, &Tmp.matrix, 0.0, &Kt.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    return cum_status;

}




ssm_err_code_t ssm_kalman_update(ssm_X_t *X, ssm_row_t *row, double t, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *like){

    int status;
    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
    gsl_vector_view pred_error = gsl_vector_subvector(&calc->_pred_error[0],0,row->ts_nonan_length);
    gsl_matrix_view Kt = gsl_matrix_submatrix(calc->_Kt,0,0,m,row->ts_nonan_length);
    gsl_matrix_view Tmp = gsl_matrix_submatrix(calc->_Tmp_N_TS_N_SV,0,0,row->ts_nonan_length,m);
    gsl_vector_view X_sv = gsl_vector_view_array(X->proj,m);
    gsl_matrix_view Ht = gsl_matrix_submatrix(calc->_Ht,0,0,m,row->ts_nonan_length);
    gsl_matrix_view Ct =  gsl_matrix_view_array(&X->proj[m], m, m);

    ssm_err_code_t cum_status = ssm_kalman_gain_computation(row, t, X, par, calc, nav);

    //////////////////
    // state update //
    //////////////////
    // X_sv += Kt * pred_error
    status = gsl_blas_dgemv(CblasNoTrans,1.0,&Kt.matrix,&pred_error.vector,1.0,&X_sv.vector);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    ///////////////////////
    // covariance update //
    ///////////////////////
    // Ct = Ct - Kt * Ht' * Ct
    status = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Ht.matrix, &Ct.matrix, 0.0, &Tmp.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &Kt.matrix, &Tmp.matrix, 0.0, &Ct.matrix);
    cum_status |=  (status != GSL_SUCCESS) ? SSM_ERR_KAL : SSM_SUCCESS;

    return cum_status;
}



/**
 * For eval_jac
 */
double ssm_diff_derivative(double jac_tpl, ssm_X_t *X, ssm_nav_t *nav, int ind)
{

    /*
      Rational basis
      we have an equation (for instance dI/dt) named eq and let's say that we are interested in its derivative against v (we assume that v follows a diffusion)'
      The template gives us d eq/d v (jac_tpl)
      However, as v can be transform (let's say log here) we want d eq / d log(v)
      The chain rule gives us:
      d eq/ dv = d eq / d log(v) * d log(v)/dv = jac_tpl
      so
      d eq / d log(v) = ( d eq / dv ) / ( d log(v) / dv)

      so in term of C:
      d eq / d log(v) = jac_tpl / r->f_derivative(v, ..)
      jac_der is the C term of v, provided by the template

      As v (jac_der) is in the scale of s_par, in case of logit_ab
      transfo, we need to provide a and b in the scale of s_par, that
      is router min_z and max_z
     */

    ssm_state_t *p = nav->states_diff->p[ind];
    
    if(jac_tpl){
        return jac_tpl / p->f_derivative(X->proj[p->offset]);
    }

    return 0.0;
}
