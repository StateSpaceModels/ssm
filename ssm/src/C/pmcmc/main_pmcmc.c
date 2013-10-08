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

int main(int argc, char *argv[])
{
    ssm_options_t *opts = ssm_options_new();
    ssm_load_options(opts, SSM_PMCMC, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav, data, fitness, opts);
    ssm_X_t ***D_J_X = ssm_D_J_X_new(data, fitness, nav, opts);
    ssm_X_t ***D_J_X_tmp = ssm_D_J_X_new(data, fitness, nav, opts);
    ssm_X_t **D_X = ssm_D_X_new(data, nav, opts); //to store sampled trajectories
    ssm_X_t **D_X_prev = ssm_D_X_new(data, nav, opts);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc[0], nav);

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
    int j, n;
    int m = 0;
    ssm_err_code_t success = SSM_SUCCESS;

    ssm_par2X(D_J_X[0][0], par, calc[0], nav);
    for(j=1; j<fitness->J; j++){
        ssm_X_copy(D_J_X[0][j], D_J_X[0][0]);
    }

    //TODO: RUN SMC success |= run_smc(...)

    if(success != SSM_SUCCESS){
        ssm_print_err("epic fail, initialization step failed");
	exit(EXIT_FAILURE);
    }

    //the first run is accepted
    fitness->log_like_prev = fitness->log_like;

    if ( ( nav->print & SSM_PRINT_X_SMOOTH ) && data->n_obs ) {
	ssm_sample_traj(D_X, D_J_X, calc[0], data, fitness);
	for(n=0; n<data->n_obs; n++){
	    ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
	    ssm_print_X(stdout, D_X_prev[n+1], par, nav, calc[0], data->rows[n], m);
	}
    }

    if(nav->print & SSM_PRINT_TRACE){
	ssm_print_trace(stdout, theta, nav, fitness->log_like_prev, m);
    }

    ////////////////
    // iterations //
    ////////////////
    double sd_fac;
    double ratio;
    for(m=1; m<n_iter; m++) {

	success = SSM_SUCCESS;

	var = ssm_adapt_eps_var_sd_fac(&sd_fac, adapt, var_input, nav, m);

	ssm_theta_ran(proposed, theta, var, sd_fac, calc[0], nav, 1);
	ssm_theta2input(input, proposed, nav);
	ssm_input2par(par, input, calc[0], nav);

	success |= ssm_check_ic(par, calc[0]);

	if(success == SSM_SUCCESS){
	    ssm_par2X(D_J_X[0][0], par, calc[0], nav);
	    D_J_X[0][0]->dt = D_J_X[0][0]->dt0;

	    for(j=1; j<fitness->J; j++){
		ssm_X_copy(D_J_X[0][j], D_J_X[0][0]);
	    }

	    //TODO: RUN SMC success |= run_smc(...)
	    success |=  ssm_metropolis_hastings(&ratio, proposed, theta, var, sd_fac, fitness, nav, calc[0], 1);
	}

	if(success == SSM_SUCCESS){ //everything went well and the proposed theta was accepted 
	    fitness->log_like_prev = fitness->log_like;
	    ssm_theta_copy(theta, proposed);

	    if ( (nav->print & SSM_PRINT_X_SMOOTH) && data->n_obs ) {
		ssm_sample_traj(D_X, D_J_X, calc[0], data, fitness);
		for(n=0; n<data->n_obs; n++){    
		    ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
		}
	    }
	} else { //failure or rejection
	    //reset par (currently corresponding to proposed) from the previous theta (theta as opposed to proposed) so that the prints (X, sample_traj) got the right values
	    ssm_theta2input(input, theta, nav);
	    ssm_input2par(par, input, calc[0], nav);
	}

	ssm_adapt_ar(adapt, is_accepted, m); //compute acceptance rate
	ssm_adapt_var(adapt, theta, m);  //compute empirical variance

	if ( (nav->print & SSM_PRINT_X_SMOOTH) && ( (m % thin_traj) == 0) ) {
	    for(n=0; n<data->n_obs; n++){
		ssm_print_X(stdout, D_X_prev[n+1], par, nav, calc[0], data->rows[n], m);
	    }
	}

	if (nav->print & SSM_PRINT_TRACE){
	    ssm_print_trace(stdout, theta, nav, fitness->log_like_prev, m);
	}

	if (nav->print & SSM_PRINT_ACC) {
	    ssm_print_ar(stdout, adapt, m);
	}

    }


    json_decref(jparameters);

    ssm_D_J_X_free(D_J_X, data, fitness);
    ssm_D_J_X_free(D_J_X_tmp, data, fitness);
    ssm_D_X_free(D_X, data);
    ssm_D_X_free(D_X_prev, data);

    ssm_N_calc_free(calc, nav);

    ssm_data_free(data);
    ssm_nav_free(nav);

    ssm_fitness_free(fitness);

    ssm_input_free(input);
    ssm_par_free(par);

    ssm_theta_free(theta);
    ssm_theta_free(proposed);
    ssm_var_free(var_input);
    ssm_adapt_free(adapt);

    return 0;
}
