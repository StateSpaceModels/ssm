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
    int i, j, n, t0, t1, id;
    char str[SSM_STR_BUFFSIZE];

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

    void *context = NULL, *sender = NULL, *receiver = NULL, *controller = NULL;
    pthread_t *workers = NULL;
    ssm_params_worker_inproc_t *params = NULL;

    if(calc[0]->threads_length >1){
        context = zmq_ctx_new();

        sender = zmq_socket (context, ZMQ_PUSH);
        snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_sender_%d", opts->id);
        zmq_bind (sender, str);

        receiver = zmq_socket (context, ZMQ_PULL);
        snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_receiver_%d", opts->id);
        zmq_bind (receiver, str);

        controller = zmq_socket (context, ZMQ_PUB);
        snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_controller_%d", opts->id);
        zmq_bind (controller, str);

        workers = malloc(calc[0]->threads_length * sizeof (pthread_t));
        params =  malloc(calc[0]->threads_length * sizeof (ssm_params_worker_inproc_t));
        int J_chunk = fitness->J / calc[0]->threads_length;
        for(i=0; i<calc[0]->threads_length; i++){
            params[i].id = opts->id;
            params[i].context = context;
	    params[i].compute_fitness = 1;
	    params[i].is_J_par = 0;
	    params[i].is_D_J_X = 0;	   
            params[i].J_chunk = J_chunk;
            params[i].data = data;
            params[i].J_par = &par;
            params[i].D_J_X = &J_X;
            params[i].calc = calc[i];
            params[i].nav = nav;
            params[i].fitness = fitness;
            params[i].f_pred = f_pred;

            pthread_create(&workers[i], NULL, ssm_worker_inproc, (void*) &params[i]);
        }

        //wait that all worker are connected
        for (i = 0; i < calc[0]->threads_length; i++) {
            zmq_recv(receiver, &id, sizeof (int), 0);
        }
    }

    for(n=0; n<data->n_obs; n++) {
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;

        if(calc[0]->threads_length > 1){

            //send work
            for (i=0; i<calc[0]->threads_length; i++) {
                zmq_send(sender, &i, sizeof (int), ZMQ_SNDMORE);
                zmq_send(sender, &n, sizeof (int), 0);
            }

            //get results from the workers
            for (i=0; i<calc[0]->threads_length; i++) {
                zmq_recv(receiver, &id, sizeof (int), 0);
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

    ssm_pipe_theta(stdout, jparameters, theta, NULL, nav);

    json_decref(jparameters);

    if(calc[0]->threads_length >1){
        zmq_send (controller, "KILL", 5, 0);
        zmq_close (sender);
        zmq_close (receiver);
        zmq_close (controller);

        for(i = 0; i < calc[0]->threads_length; i++){
            pthread_join(workers[i], NULL);
        }

        free(workers);
        free(params);
        zmq_ctx_destroy (context);
    }

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
