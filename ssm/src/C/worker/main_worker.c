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
    int j, n, t0, t1;
    char str[SSM_STR_BUFFSIZE];

    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_WORKER, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);
    opts->J = 1;

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t *calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);
    ssm_X_t *X = ssm_X_new(nav, opts);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc, nav);

    void *context = zmq_ctx_new();

    // Socket to server controller
    void *server_controller = zmq_socket (context, ZMQ_SUB);
    snprintf(str, SSM_STR_BUFFSIZE, "tcp://%s:%d", opts->server, 5559);
    zmq_connect (server_controller, str);
    zmq_setsockopt (server_controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *server_receiver = zmq_socket (context, ZMQ_PULL);
    snprintf(str, SSM_STR_BUFFSIZE, "tcp://%s:%d", opts->server, 5557);
    zmq_connect (server_receiver, str);

    //  Socket to send messages (results) to the server
    void *server_sender = zmq_socket (context, ZMQ_PUSH);
    snprintf(str, SSM_STR_BUFFSIZE, "tcp://%s:%d", opts->server, 5558);
    zmq_connect (server_sender, str);
   
    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    zmq_pollitem_t items [] = {
        { server_receiver, 0, ZMQ_POLLIN, 0 },
        { server_controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {
	    
            //get a particle from the server
	    zmq_recv(server_receiver, &n, sizeof (int), 0);
	    ssm_zmq_recv_par(par, server_receiver);
	    zmq_recv(server_receiver, &j, sizeof (int), 0);
	    //printf("j: %d j %d\n", j, j);
	    ssm_zmq_recv_X(X, server_receiver);
	    zmq_recv(server_receiver, &(fitness->cum_status[0]), sizeof (ssm_err_code_t), 0);

	    //do the computations..
            t0 = (n) ? data->rows[n-1]->time: 0;
            t1 = data->rows[n]->time;
	    ssm_X_reset_inc(X, data->rows[n], nav);
	    fitness->cum_status[0] |= (*f_pred)(X, t0, t1, par, nav, calc);
	    if((opts->worker_algo != SSM_SIMUL) && data->rows[n]->ts_nonan_length) {
		fitness->weights[0] = (fitness->cum_status[0] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], X, par, calc, nav, fitness)) : 0.0;
		fitness->cum_status[0] = SSM_SUCCESS;
	    }

	    //send results
	    zmq_send(server_sender, &j, sizeof (int), ZMQ_SNDMORE);    
	    ssm_zmq_send_X(server_sender, X, ZMQ_SNDMORE);
	    if(opts->worker_algo != SSM_SIMUL){
		zmq_send(server_sender, &(fitness->weights[0]), sizeof (double), ZMQ_SNDMORE);
	    }
	    zmq_send(server_sender, &(fitness->cum_status[0]), sizeof (ssm_err_code_t), 0);
	    //printf("j: %d j %d sent back\n", j, j);

        }

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf[SSM_STR_BUFFSIZE];
	    zmq_recv(server_controller, buf, SSM_STR_BUFFSIZE, 0);           

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (server_receiver);
    zmq_close (server_sender);
    zmq_close (server_controller);

    ssm_X_free(X);
    ssm_calc_free(calc, nav);
    ssm_data_free(data);
    ssm_fitness_free(fitness);

    ssm_nav_free(nav);
    ssm_input_free(input);
    ssm_par_free(par);

    zmq_ctx_destroy(context);

    return 0;
}
