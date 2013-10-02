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

/**
   used for euler multinomial integrarion. When duration of
   infection is close to the time step duration, the method becomes
   inacurate (the waiting time is geometric instead of
   exponential. So we ensure that the rate has the correct magnitude
   by correcting it
*/
double correct_rate(double rate, double dt)
{
    return -log(1.0-rate*dt)/dt;
}


/**
 * Check if remainder has not become negative
 */
ssm_err_code_t ssm_check_no_neg_remainder(ssm_X_t *p_X, ssm_nav_t *nav, ssm_calc_t *calc, double t)
{
    int i;
    ssm_it_states_t *rem = nav->states_remainders;

    for(i=0; i<rem->length; i++){
        if (rem->p[i]->f_remainder(p_X, calc, t) < 0.0){
            return SSM_ERR_REM;
        }
    }

    return SSM_SUCCESS;
}


ssm_f_pred_t ssm_get_f_pred(ssm_calc_t *calc)
{
    ssm_implementations_t implementation = calc->implementation;
    ssm_noises_off_t noises_off= calc->noises_off;

    if (implementation == PLOM_ODE) {
        return &ssm_f_prediction_ode;

    } else if (implementation == PLOM_SDE){

        if (noises_off == (SSM_NO_DEM_STO | SSM_NO_WHITE_NOISE | SSM_NO_DIFF) ) {
            return &ssm_f_prediction_ode;
        } else if (noises_off == (SSM_NO_DEM_STO | SSM_NO_WHITE_NOISE) ) {
            return &ssm_f_prediction_sde_no_dem_sto_no_white_noise;
        } else if (noises_off == (SSM_NO_DEM_STO | SSM_NO_DIFF) ) {
            return &ssm_f_prediction_sde_no_dem_sto_no_diff;
        } else if (noises_off == (SSM_NO_WHITE_NOISE | SSM_NO_DIFF) ) {
            return &ssm_f_prediction_sde_no_white_noise_no_diff;
        } else if (noises_off == SSM_NO_DEM_STO ) {
            return &ssm_f_prediction_sde_no_dem_sto;
        } else if (noises_off == SSM_NO_WHITE_NOISE ) {
            return &ssm_f_prediction_sde_no_white_noise;
        } else if (noises_off == SSM_NO_DIFF ) {
            return &ssm_f_prediction_sde_no_diff;
        } else {
            return &ssm_f_prediction_sde_full;
        }

    } else if (implementation == PLOM_PSR){
        //no_sto_env is handled within the step funciton
        if(noises_off &ssm_ SSM_NO_DIFF){
            return &ssm_f_prediction_psr_no_diff;
        } else {
            return &ssm_f_prediction_psr;
        }
    }

    return NULL;
}


ssm_err_code_t ssm_f_prediction_ode(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t=t0;
    double h = p_X->dt; //h is the initial integration step size
    calc->_par = par; //pass the ref to par so that it is available wihtin the function to integrate

    double *y = p_X->proj;

    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply (calc->evolve, calc->control, calc->step, &(calc->sys), &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            return SSM_ERR_ODE;
        }
    }

    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}



ssm_err_code_t ssm_f_prediction_sde_no_dem_sto_no_white_noise(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_dem_sto_no_white_noise(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }

    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_sde_no_dem_sto_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_dem_sto(p_X, t, par, nav, calc);
        t += p_X->dt;
    }

    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_sde_no_white_noise_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_white_noise(p_X, t, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_sde_no_dem_sto(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_dem_sto(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_sde_no_white_noise(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_white_noise(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_sde_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_full(p_X, t, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_sde_full(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_full(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_ekf_no_dem_sto_no_white_noise(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t=t0;
    double h = p_X->dt; //h is the initial integration step size
    calc->_par = par; //pass the ref to par so that it is available wihtin the function to integrate

    double *y = p_X->proj;

    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply (calc->evolve, calc->control, calc->step, &(calc->sys), &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            return SSM_ERR_ODE;
        }
    }

    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_ekf_no_dem_sto_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_dem_sto(p_X, t, par, nav, calc);
        t += p_X->dt;
    }

    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_ekf_no_white_noise_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_white_noise(p_X, t, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_ekf_no_dem_sto(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_dem_sto(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_ekf_no_white_noise(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_no_white_noise(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_ekf_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_full(p_X, t, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}

ssm_err_code_t ssm_f_prediction_ekf_full(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_sde_full(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


ssm_err_code_t ssm_f_prediction_psr(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_psr(p_X, t, par, nav, calc);
        ssm_compute_diff(p_X, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}



ssm_err_code_t ssm_f_prediction_psr_no_diff(ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    double t = t0;

    while (t < t1) {
        ssm_step_psr(p_X, t, par, nav, calc);
        t += p_X->dt;
    }
    return ssm_check_no_neg_remainder(p_X, nav, calc, t1);
}


