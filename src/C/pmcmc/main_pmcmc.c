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

static ssm_err_code_t run_smc(ssm_err_code_t (*f_pred) (ssm_X_t *, double, double, ssm_par_t *, ssm_nav_t *, ssm_calc_t *), ssm_X_t ***D_J_X, ssm_X_t ***D_J_X_tmp, ssm_par_t *par, ssm_calc_t **calc, ssm_data_t *data, ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_workers_t *workers)
{
    int i, j, n, np1, id, the_j;
    double t0, t1;

    fitness->log_like = 0.0;
    fitness->log_prior = 0.0;
    fitness->n_all_fail = 0;

    for(j=0; j<fitness->J; j++){
	fitness->cum_status[j] = SSM_SUCCESS;
    }

    for(n=0; n<data->n_obs; n++) {
        np1 = n+1;
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;

	if(!workers->flag_tcp){
	    for(j=0; j<fitness->J; j++){
		ssm_X_copy(D_J_X[np1][j], D_J_X[n][j]);
	    }
	}

	if(workers->flag_tcp){
	    //send work
	    for (j=0;j<fitness->J;j++) {
		zmq_send(workers->sender, &n, sizeof (int), ZMQ_SNDMORE);
		ssm_zmq_send_par(workers->sender, par, ZMQ_SNDMORE);

		zmq_send(workers->sender, &j, sizeof (int), ZMQ_SNDMORE);                   	       	       
		ssm_zmq_send_X(workers->sender, D_J_X[n][j], ZMQ_SNDMORE);
		zmq_send(workers->sender, &(fitness->cum_status[j]), sizeof (ssm_err_code_t), 0);
	    }

	    //get results from the workers
	    for (j=0; j<fitness->J; j++) {
		zmq_recv(workers->receiver, &the_j, sizeof (int), 0);
		ssm_zmq_recv_X(D_J_X[ np1 ][ the_j ], workers->receiver);
		zmq_recv(workers->receiver, &(fitness->weights[the_j]), sizeof (double), 0);
		zmq_recv(workers->receiver, &(fitness->cum_status[the_j]), sizeof (ssm_err_code_t), 0);
	    }

	} else if(calc[0]->threads_length > 1){
	    //send work
            for (i=0; i<calc[0]->threads_length; i++) {
                zmq_send(workers->sender, &i, sizeof (int), ZMQ_SNDMORE);
                zmq_send(workers->sender, &n, sizeof (int), 0);
            }

            //get results from the workers
            for (i=0; i<calc[0]->threads_length; i++) {
                zmq_recv(workers->receiver, &id, sizeof (int), 0);
            }
        } else {

	    for(j=0;j<fitness->J;j++) {
		ssm_X_reset_inc(D_J_X[np1][j], data->rows[n], nav);
		fitness->cum_status[j] |= (*f_pred)(D_J_X[np1][j], t0, t1, par, nav, calc[0]);
		if(data->rows[n]->ts_nonan_length) {
		    fitness->weights[j] = (fitness->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], D_J_X[np1][j], par, calc[0], nav, fitness)) : 0.0;
		    fitness->cum_status[j] = SSM_SUCCESS;
		}
	    }
	}
	
        if(data->rows[n]->ts_nonan_length) {
            if(ssm_weight(fitness, data->rows[n], nav, n)) {
                ssm_systematic_sampling(fitness, calc[0], n);
            }
            ssm_resample_X(fitness, &D_J_X[np1], &D_J_X_tmp[np1], n);
        }
    }
    return ( (data->n_obs != 0) && (fitness->n_all_fail == data->n_obs) ) ? SSM_ERR_PRED: SSM_SUCCESS;
}

int main(int argc, char *argv[])
{
    char str[SSM_STR_BUFFSIZE];

    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_PMCMC, argc, argv);

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
    ssm_par_t *par_proposed = ssm_par_new(input, calc[0], nav);

    ssm_theta_t *theta = ssm_theta_new(input, nav);
    ssm_theta_t *proposed = ssm_theta_new(input, nav);
    ssm_var_t *var_input = ssm_var_new(jparameters, nav);
    ssm_var_t *var; //the covariance matrix used;
    ssm_adapt_t *adapt = ssm_adapt_new(nav, opts);

    int n_iter = opts->n_iter;
    int n_traj = GSL_MIN(n_iter, opts->n_traj);
    int thin_traj = (int) ( (double) n_iter / (double) n_traj); //the thinning interval

    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    ssm_workers_t *workers = ssm_workers_start(D_J_X, &par, data, calc, fitness, f_pred, nav, opts, SSM_WORKER_D_X | SSM_WORKER_FITNESS);

    /////////////////////////
    // initialization step //
    /////////////////////////
    int j, n;
    int m = 0;

    ssm_par2X(D_J_X[0][0], par, calc[0], nav);
    for(j=1; j<fitness->J; j++){
        ssm_X_copy(D_J_X[0][j], D_J_X[0][0]);
    }

    ssm_err_code_t success = run_smc(f_pred, D_J_X, D_J_X_tmp, par_proposed, calc, data, fitness, nav, workers);
    success |= ssm_log_prob_prior(&fitness->log_prior, proposed, nav, fitness);

    if(success != SSM_SUCCESS){
        ssm_print_err("epic fail, initialization step failed");
        exit(EXIT_FAILURE);
    }

    //the first run is accepted
    fitness->log_like_prev = fitness->log_like;
    fitness->log_prior_prev = fitness->log_prior;

    if ( ( nav->print & SSM_PRINT_X ) && data->n_obs ) {
        ssm_sample_traj(D_X, D_J_X, calc[0], data, fitness);
        for(n=0; n<data->n_obs; n++){
            ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
            ssm_print_X(nav->X, D_X_prev[n+1], par, nav, calc[0], data->rows[n], m);
        }
    }

    if(nav->print & SSM_PRINT_TRACE){
        ssm_print_trace(nav->trace, theta, nav, fitness->log_like_prev + fitness->log_prior_prev, m);
    }

    ssm_dic_init(fitness, fitness->log_like_prev);

    if (nav->print & SSM_PRINT_LOG) {
	snprintf(str, SSM_STR_BUFFSIZE, "%d\t logLike.: %g\t accepted: %d\t acc. rate: %g", m, fitness->log_like_prev + fitness->log_prior_prev, !(success & SSM_MH_REJECT), adapt->ar);
	ssm_print_log(str);
    }

    ////////////////
    // iterations //
    ////////////////
    double sd_fac;
    double ratio;
    for(m=1; m<n_iter; m++) {
        var = ssm_adapt_eps_var_sd_fac(&sd_fac, adapt, var_input, nav, m);
        ssm_theta_ran(proposed, theta, var, sd_fac, calc[0], nav, 1);
        ssm_theta2input(input, proposed, nav);
        ssm_input2par(par_proposed, input, calc[0], nav);

        success = ssm_check_ic(par_proposed, calc[0]);

        if(success == SSM_SUCCESS){
            ssm_par2X(D_J_X[0][0], par_proposed, calc[0], nav);
            D_J_X[0][0]->dt = D_J_X[0][0]->dt0;
            for(j=1; j<fitness->J; j++){
                ssm_X_copy(D_J_X[0][j], D_J_X[0][0]);
            }

	    success |= run_smc(f_pred, D_J_X, D_J_X_tmp, par_proposed, calc, data, fitness, nav, workers);
            success |= ssm_metropolis_hastings(fitness, &ratio, proposed, theta, var, sd_fac, nav, calc[0], 1);
        }

        if(success == SSM_SUCCESS){ //everything went well and the proposed theta was accepted
            fitness->log_like_prev = fitness->log_like;
            fitness->log_prior_prev = fitness->log_prior;
            ssm_theta_copy(theta, proposed);
            ssm_par_copy(par, par_proposed);

            if ( (nav->print & SSM_PRINT_X) && data->n_obs ) {
                ssm_sample_traj(D_X, D_J_X, calc[0], data, fitness);
                for(n=0; n<data->n_obs; n++){
                    ssm_X_copy(D_X_prev[n+1], D_X[n+1]);
                }
            }
        }

        ssm_adapt_ar(adapt, (success == SSM_SUCCESS) ? 1: 0, m); //compute acceptance rate
        ssm_adapt_var(adapt, theta, m);  //compute empirical variance

        if ( (nav->print & SSM_PRINT_X) && ( (m % thin_traj) == 0) ) {
            for(n=0; n<data->n_obs; n++){
                ssm_print_X(nav->X, D_X_prev[n+1], par, nav, calc[0], data->rows[n], m);
            }
        }

        if (nav->print & SSM_PRINT_TRACE){
            ssm_print_trace(nav->trace, theta, nav, fitness->log_like_prev + fitness->log_prior_prev, m);
        }
	ssm_dic_update(fitness, fitness->log_like_prev);


        if (nav->print & SSM_PRINT_DIAG) {
            ssm_print_ar(nav->diag, adapt, m);
        }

	if (nav->print & SSM_PRINT_LOG) {
	    snprintf(str, SSM_STR_BUFFSIZE, "%d\t logLike.: %g\t accepted: %d\t acc. rate: %g", m, fitness->log_like_prev + fitness->log_prior_prev, !(success & SSM_MH_REJECT), adapt->ar);
	    ssm_print_log(str);
	}
    }

    if (!(nav->print & SSM_PRINT_LOG)) {
	ssm_dic_end(fitness, nav, m);
	ssm_pipe_theta(stdout, jparameters, theta, var, fitness, nav, opts);
    }

    json_decref(jparameters);

    ssm_workers_stop(workers);

    ssm_D_J_X_free(D_J_X, data, fitness);
    ssm_D_J_X_free(D_J_X_tmp, data, fitness);
    ssm_D_X_free(D_X, data);
    ssm_D_X_free(D_X_prev, data);

    ssm_N_calc_free(calc, nav);

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
