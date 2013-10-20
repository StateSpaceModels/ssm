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

static ssm_err_code_t run_kalman_and_store_traj(ssm_X_t **D_X, ssm_par_t *par, ssm_fitness_t *fitness, ssm_data_t *data, ssm_calc_t *calc, ssm_nav_t *nav)
{
    int n, np1;
    double t0, t1;
    ssm_err_code_t rc;

    fitness->log_like = 0.0;
    fitness->log_prior = 0.0;

    for(n=0; n<data->n_obs; n++) {
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;
        np1 = n+1;

        ssm_X_copy(D_X[np1], D_X[n]);
        ssm_X_reset_inc(D_X[np1], data->rows[n], nav);

        rc = ssm_f_prediction_ode(D_X[np1], t0, t1, par, nav, calc);
        if(rc != SSM_SUCCESS){
            return rc;
        }
        if(data->rows[n]->ts_nonan_length) {
            rc = ssm_kalman_update(fitness, D_X[np1], data->rows[n], t1, par, calc, nav);
            if(rc != SSM_SUCCESS){
                return rc;
            }
        }
    }

    return SSM_SUCCESS;
}


int main(int argc, char *argv[])
{
    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_KMCMC, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t *calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);
    ssm_X_t **D_X = ssm_D_X_new(data, nav, opts); //to store trajectory
    ssm_X_t **D_X_prev = ssm_D_X_new(data, nav, opts);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc, nav);
    ssm_par_t *par_proposed = ssm_par_new(input, calc, nav);

    ssm_theta_t *theta = ssm_theta_new(input, nav);
    ssm_theta_t *proposed = ssm_theta_new(input, nav);
    ssm_var_t *var_input = ssm_var_new(jparameters, nav);
    ssm_var_t *var; //the covariance matrix used;
    ssm_adapt_t *adapt = ssm_adapt_new(nav, opts);

    int n_iter = opts->n_iter;
    int n_traj = GSL_MIN(n_iter, opts->n_traj);
    int thin_traj = (int) ( (double) n_iter / (double) n_traj); //the thinning interval


    /////////////////////////
    // initialization step //
    /////////////////////////
    int n;
    int m = 0;

    ssm_par2X(D_X[0], par, calc, nav);

    ssm_err_code_t success = run_kalman_and_store_traj(D_X, par, fitness, data, calc, nav);
    success |= ssm_log_prob_prior(&fitness->log_prior, proposed, nav, fitness);
    if(success != SSM_SUCCESS){
        ssm_print_err("epic fail, initialization step failed");
        exit(EXIT_FAILURE);
    }

    //the first run is accepted
    fitness->log_like_prev = fitness->log_like;
    fitness->log_prior_prev = fitness->log_prior;

    if ( ( nav->print & SSM_PRINT_X ) && data->n_obs ) {
        for(n=0; n<data->n_obs; n++){
            ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
            ssm_print_X(nav->X, D_X_prev[n+1], par, nav, calc, data->rows[n], m);
        }
    }

    if(nav->print & SSM_PRINT_TRACE){
        ssm_print_trace(nav->trace, theta, nav, fitness->log_like_prev + fitness->log_prior_prev, m);
    }

    ////////////////
    // iterations //
    ////////////////
    double sd_fac;
    double ratio;
    for(m=1; m<n_iter; m++) {
        success = SSM_SUCCESS;

        var = ssm_adapt_eps_var_sd_fac(&sd_fac, adapt, var_input, nav, m);

        ssm_theta_ran(proposed, theta, var, sd_fac, calc, nav, 1);
        ssm_theta2input(input, proposed, nav);
        ssm_input2par(par_proposed, input, calc, nav);

        success |= ssm_check_ic(par_proposed, calc);

        if(success == SSM_SUCCESS){
            ssm_par2X(D_X[0], par_proposed, calc, nav);
            D_X[0]->dt = D_X[0]->dt0;
            ssm_kalman_reset_Ct(D_X[0], nav);

            success |= run_kalman_and_store_traj(D_X, par, fitness, data, calc, nav);
            success |= ssm_metropolis_hastings(fitness, &ratio, proposed, theta, var, sd_fac, nav, calc, 1);
        }

        if(success == SSM_SUCCESS){ //everything went well and the proposed theta was accepted
            fitness->log_like_prev = fitness->log_like;
            fitness->log_prior_prev = fitness->log_prior;
            ssm_theta_copy(theta, proposed);
            ssm_par_copy(par, par_proposed);

            if ( (nav->print & SSM_PRINT_X) && data->n_obs ) {
                for(n=0; n<data->n_obs; n++){
                    ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
                }
            }
        }

        ssm_adapt_ar(adapt, (success == SSM_SUCCESS) ? 1: 0, m); //compute acceptance rate
        ssm_adapt_var(adapt, theta, m);  //compute empirical variance

        if ( (nav->print & SSM_PRINT_X) && ( (m % thin_traj) == 0) ) {
            for(n=0; n<data->n_obs; n++){
                ssm_print_X(nav->X, D_X_prev[n+1], par, nav, calc, data->rows[n], m);
            }
        }

        if (nav->print & SSM_PRINT_TRACE){
            ssm_print_trace(nav->trace, theta, nav, fitness->log_like_prev + fitness->log_prior_prev, m);
        }

        if (nav->print & SSM_PRINT_DIAG) {
            ssm_print_ar(nav->diag, adapt, m);
        }
    }

    if (!(nav->print & SSM_PRINT_LOG)) {
	ssm_pipe_theta(stdout, jparameters, theta, var, nav, opts);
    }

    json_decref(jparameters);

    ssm_D_X_free(D_X, data);
    ssm_D_X_free(D_X_prev, data);

    ssm_calc_free(calc, nav);

    ssm_data_free(data);
    ssm_nav_free(nav);

    ssm_fitness_free(fitness);

    ssm_input_free(input);
    ssm_par_free(par_proposed);
    ssm_par_free(par);

    ssm_theta_free(theta);
    ssm_theta_free(proposed);
    ssm_var_free(var_input);
    ssm_adapt_free(adapt);

    return 0;
}
