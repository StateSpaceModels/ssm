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
    int i, j, n, t0, t1, id, the_j;

    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_SMC, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);

    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav, data, fitness, opts);
    ssm_X_t **J_X = ssm_J_X_new(fitness, nav, opts);
    ssm_X_t **J_X_tmp = ssm_J_X_new(fitness, nav, opts);
    ssm_hat_t *hat = ssm_hat_new(nav);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc[0], nav);
    ssm_theta_t *theta = ssm_theta_new(input, nav);

    int flag_prior = opts->flag_prior;
    int flag_no_filter = opts->flag_no_filter;

    ssm_par2X(J_X[0], par, calc[0], nav);
    for(j=1; j<fitness->J; j++){
        ssm_X_copy(J_X[j], J_X[0]);
    }

    for(j=0; j<fitness->J; j++) {
        fitness->cum_status[j] = SSM_SUCCESS;
    }

    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    ssm_workers_t *workers = ssm_workers_start(&J_X, &par, data, calc, fitness, f_pred, nav, opts, SSM_WORKER_FITNESS);

    for(n=0; n<data->n_obs; n++) {
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;

	if(workers->flag_tcp){
	    //send work
	    for (j=0;j<fitness->J;j++) {
		zmq_send(workers->sender, &n, sizeof (int), ZMQ_SNDMORE);
		ssm_zmq_send_par(workers->sender, par, ZMQ_SNDMORE);

		zmq_send(workers->sender, &j, sizeof (int), ZMQ_SNDMORE);                   	       	       
		ssm_zmq_send_X(workers->sender, J_X[j], ZMQ_SNDMORE);
		zmq_send(workers->sender, &(fitness->cum_status[j]), sizeof (ssm_err_code_t), 0);
		//printf("part %d sent %d\n", j, 0);
	    }

	    //get results from the workers
	    for (j=0; j<fitness->J; j++) {
		zmq_recv(workers->receiver, &the_j, sizeof (int), 0);
		ssm_zmq_recv_X(J_X[ the_j ], workers->receiver);
		zmq_recv(workers->receiver, &(fitness->weights[the_j]), sizeof (double), 0);
		zmq_recv(workers->receiver, &(fitness->cum_status[the_j]), sizeof (ssm_err_code_t), 0);
		//printf("part  %d received\n", the_j);
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
                ssm_X_reset_inc(J_X[j], data->rows[n], nav);
                fitness->cum_status[j] |= (*f_pred)(J_X[j], t0, t1, par, nav, calc[0]);

                if(data->rows[n]->ts_nonan_length) {
                    fitness->weights[j] = (fitness->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], J_X[j], par, calc[0], nav, fitness)) : 0.0;
                    fitness->cum_status[j] = SSM_SUCCESS;
                }
            }

        }

        if(!flag_no_filter && data->rows[n]->ts_nonan_length) {
            if(ssm_weight(fitness, data->rows[n], nav, n)) {
		ssm_systematic_sampling(fitness, calc[0], n);
            }

            if (nav->print & SSM_PRINT_HAT) {
                ssm_hat_eval(hat, J_X, &par, nav, calc[0], fitness, t1, 0);
            }

            if (nav->print & SSM_PRINT_DIAG) {
                ssm_print_pred_res(nav->diag, J_X, par, nav, calc[0], data, data->rows[n], fitness);
            }

	    ssm_resample_X(fitness, &J_X, &J_X_tmp, n);

        } else if (nav->print & SSM_PRINT_HAT) { //we do not filter or all data ara NaN (no info).
            ssm_hat_eval(hat, J_X, &par, nav, calc[0], NULL, t1, 0);
        }

        if (nav->print & SSM_PRINT_HAT) {
            ssm_print_hat(nav->hat, hat, nav, data->rows[n]);
        }

        if (nav->print & SSM_PRINT_X) {
            for(j=0; j<fitness->J; j++) {
                ssm_print_X(nav->X, J_X[j], par, nav, calc[0], data->rows[n], j);
            }
        }
    }

    if (flag_prior) {
        double log_prob_prior_value;
        ssm_err_code_t rc = ssm_log_prob_prior(&log_prob_prior_value, theta, nav, fitness);
        if(rc != SSM_SUCCESS && (nav->print & SSM_PRINT_WARNING)){
            ssm_print_warning("error log_prob_prior computation");
        }
        fitness->log_like += log_prob_prior_value;
    }

    if (nav->print & SSM_PRINT_TRACE) {
        ssm_print_trace(nav->trace, theta, nav, fitness->log_like, 0);
    }

    if (!(nav->print & SSM_PRINT_LOG)) {
	ssm_pipe_theta(stdout, jparameters, theta, NULL, nav, opts);
    } else {
	char str[SSM_STR_BUFFSIZE];
	snprintf(str, SSM_STR_BUFFSIZE, "logLike.: %g", fitness->log_like);
	ssm_print_log(str);
    }

    json_decref(jparameters);

    ssm_workers_stop(workers);

    ssm_J_X_free(J_X, fitness);
    ssm_J_X_free(J_X_tmp, fitness);
    ssm_hat_free(hat);
    ssm_N_calc_free(calc, nav);
    ssm_data_free(data);
    ssm_fitness_free(fitness);

    ssm_nav_free(nav);
    ssm_theta_free(theta);
    ssm_input_free(input);
    ssm_par_free(par);




    return 0;
}
