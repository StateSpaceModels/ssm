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
    ssm_options_load(opts, SSM_SIMUL, argc, argv);
    
    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    json_t *jresource = json_object_get(jparameters, "jresource");
    json_t *jvalues = NULL;
    for(i=0; i< json_array_size(jresource); i++){
        json_t *el = json_array_get(jresource, i);
        const char* name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "values") == 0) {
            jvalues = json_object_get(el, "data");
            break;
        }
    }

    int is_predict_from_traces = json_is_array(jvalues);
    if(is_predict_from_traces){
	opts->J = json_array_size(jvalues);
    }

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav, data, fitness, opts);
    ssm_X_t **J_X = ssm_J_X_new(fitness, nav, opts);
    ssm_hat_t *hat = ssm_hat_new(nav);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new((is_predict_from_traces) ? NULL: jparameters, nav);
    ssm_par_t **J_par = malloc(fitness->J * sizeof (ssm_par_t *));
    if(J_par == NULL) {
        ssm_print_err("Allocation impossible for ssm_par_t *");
        exit(EXIT_FAILURE);
    }

    for(j=0; j<fitness->J; j++) {
	if(is_predict_from_traces){
	    ssm_jforced(input, json_array_get(jvalues, j), nav);
	}
	J_par[j] = ssm_par_new(input, calc[0], nav);
	ssm_par2X(J_X[j], J_par[j], calc[0], nav);
    }
    
    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    ssm_workers_t *workers = ssm_workers_start(&J_X, J_par, data, calc, fitness, f_pred, nav, opts, SSM_WORKER_J_PAR);
    
    for(j=0; j<fitness->J; j++) {
        fitness->cum_status[j] = SSM_SUCCESS;
    }

    for(n=0; n<data->length; n++) {
	t0 = (n) ? data->rows[n-1]->time: 0;
	t1 = data->rows[n]->time;


	if(workers->flag_tcp){
	    //send work
	    for (j=0;j<fitness->J;j++) {
		zmq_send(workers->sender, &n, sizeof (int), ZMQ_SNDMORE);
		ssm_zmq_send_par(workers->sender, J_par[j], ZMQ_SNDMORE);

		zmq_send(workers->sender, &j, sizeof (int), ZMQ_SNDMORE);                   	       	       
		ssm_zmq_send_X(workers->sender, J_X[j], ZMQ_SNDMORE);
		zmq_send(workers->sender, &(fitness->cum_status[j]), sizeof (ssm_err_code_t), 0);
		//printf("part %d sent %d\n", j, 0);
	    }

	    //get results from the workers
	    for (j=0; j<fitness->J; j++) {
		zmq_recv(workers->receiver, &the_j, sizeof (int), 0);
		ssm_zmq_recv_X(J_X[ the_j ], workers->receiver);
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
		fitness->cum_status[j] |= f_pred(J_X[j], t0, t1, J_par[j], nav, calc[0]);
	    }

	}

	if (nav->print & SSM_PRINT_HAT) {
	    ssm_hat_eval(hat, J_X, J_par, nav, calc[0], NULL, t1, 0);
	    ssm_print_hat(nav->hat, hat, nav, data->rows[n]);
        }

	if (nav->print & SSM_PRINT_X) {
	    for(j=0; j<fitness->J; j++) {
		ssm_print_X(nav->X, J_X[j], J_par[j], nav, calc[0], data->rows[n], j);
	    }
	}
    }

    if(!(nav->print & SSM_PRINT_LOG) &&!is_predict_from_traces){
	if (!(nav->print & SSM_PRINT_HAT)) { //hat was not computed
	    ssm_hat_eval(hat, J_X, J_par, nav, calc[0], NULL, t1, 0);	
	}
	ssm_pipe_hat(stdout, jparameters, input, hat, J_par[0], calc[0], nav, opts, t1);
    }

    json_decref(jparameters);

    ssm_workers_stop(workers);

    for(j=0; j<fitness->J; j++) {
        ssm_par_free(J_par[j]);
    }
    free(J_par);

    ssm_J_X_free(J_X, fitness);
    ssm_N_calc_free(calc, nav);
    ssm_data_free(data);
    ssm_nav_free(nav);
    ssm_fitness_free(fitness);
    ssm_hat_free(hat);

    ssm_input_free(input);

    return 0;
}
