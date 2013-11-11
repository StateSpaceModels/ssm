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

void *ssm_worker_inproc(void *params)
{
    char str[SSM_STR_BUFFSIZE];
    ssm_params_worker_inproc_t *p = (ssm_params_worker_inproc_t *) params;

    int id = p->id;
    void *context = p->context;
    ssm_worker_opt_t wopts = p->wopts;
    int J_chunk = p->J_chunk;   
    ssm_data_t *data = p->data;
    ssm_par_t **J_par = p->J_par;
    ssm_X_t ***D_J_X = p->D_J_X;
    ssm_calc_t *calc = p->calc;
    ssm_nav_t *nav = p->nav;
    ssm_fitness_t *fitness = p->fitness;
    ssm_f_pred_t f_pred = p->f_pred;

    // Socket to server controller
    void *controller = zmq_socket (context, ZMQ_SUB);
    snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_controller_%d", id);
    zmq_connect (controller, str);
    zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *receiver = zmq_socket (context, ZMQ_PULL);
    snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_sender_%d", id);
    zmq_connect (receiver, str);

    //  Socket to send messages (results) to the server
    void *sender = zmq_socket (context, ZMQ_PUSH);
    snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_receiver_%d", id);
    zmq_connect (sender, str);

    // ready !
    zmq_send(sender, &calc->thread_id, sizeof (int), 0);

    zmq_pollitem_t items [] = {
        { receiver, 0, ZMQ_POLLIN, 0 },
        { controller, 0, ZMQ_POLLIN, 0 }
    };

    int j, n, t0, t1;
    int the_id;

    int _zero = 0;
    int *j_par = (SSM_WORKER_J_PAR & wopts) ? &j: &_zero;
    int np1;
    int *n_X = (SSM_WORKER_D_X & wopts) ? &np1 : &_zero;

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

            zmq_recv(receiver, &the_id, sizeof (int), 0);
            zmq_recv(receiver, &n, sizeof (int), 0);

	    np1 = n + 1;
            t0 = (n) ? data->rows[n-1]->time: 0;
            t1 = data->rows[n]->time;

            int J_start = the_id * J_chunk;
            int J_end = (the_id+1 == calc->threads_length) ? fitness->J : (the_id+1)*J_chunk;

            for(j=J_start; j<J_end; j++ ){
                ssm_X_reset_inc(D_J_X[*n_X][j], data->rows[n], nav);
		fitness->cum_status[j] |= (*f_pred)(D_J_X[*n_X][j], t0, t1, J_par[*j_par], nav, calc);
		
                if((SSM_WORKER_FITNESS & wopts) && data->rows[n]->ts_nonan_length) {
                    fitness->weights[j] = (fitness->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], D_J_X[*n_X][j], J_par[*j_par], calc, nav, fitness)) : 0.0;
                    fitness->cum_status[j] = SSM_SUCCESS;
                }
            }

            //send back id of the batch of particles now integrated
            zmq_send(sender, &calc->thread_id, sizeof (int), 0);
        }

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
            char buf[SSM_STR_BUFFSIZE];
            zmq_recv(controller, buf, SSM_STR_BUFFSIZE, 0);

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (receiver);
    zmq_close (sender);
    zmq_close (controller);

    return NULL;
}


ssm_workers_t *ssm_workers_start(ssm_X_t ***D_J_X, ssm_par_t **J_par, ssm_data_t *data, ssm_calc_t **calc, ssm_fitness_t *fitness, ssm_f_pred_t f_pred, ssm_nav_t *nav, ssm_options_t *opts, ssm_worker_opt_t wopts)
{
    int i, id;
    char str[SSM_STR_BUFFSIZE];
    
    ssm_workers_t *w = malloc(sizeof(ssm_workers_t));
    if(w == NULL){
	ssm_print_err("allocation impossible for ssm_workers_t");
	exit(EXIT_FAILURE);
    }

    w->flag_tcp = opts->flag_tcp;
    w->inproc_length = calc[0]->threads_length;
    w->wopts = wopts;

    if(opts->flag_tcp){
	w->context = zmq_ctx_new();;

        //  Socket to send messages on
        w->sender = zmq_socket(w->context, ZMQ_PUSH);
        zmq_bind(w->sender, "tcp://*:5557");

        //  Socket to receive messages on
        w->receiver = zmq_socket(w->context, ZMQ_PULL);
        zmq_bind(w->receiver, "tcp://*:5558");

        //  Socket for worker control
        w->controller = zmq_socket(w->context, ZMQ_PUB);
        zmq_bind(w->controller, "tcp://*:5559");

	w->params = NULL;
	w->workers = NULL;

    } else if (w->inproc_length == 1){
	w->context = NULL;
	w->sender = NULL;
	w->receiver = NULL;
	w->controller = NULL;
	w->params = NULL;
	w->workers = NULL;       

    } else {
	w->context = zmq_ctx_new();

	w->sender = zmq_socket (w->context, ZMQ_PUSH);
	snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_sender_%d", opts->id);
	zmq_bind (w->sender, str);

	w->receiver = zmq_socket (w->context, ZMQ_PULL);
	snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_receiver_%d", opts->id);
	zmq_bind (w->receiver, str);

	w->controller = zmq_socket (w->context, ZMQ_PUB);
	snprintf(str, SSM_STR_BUFFSIZE, "inproc://ssm_server_controller_%d", opts->id);
	zmq_bind (w->controller, str);

	w->workers = malloc(w->inproc_length * sizeof (pthread_t));
	if(w->workers == NULL){
	    ssm_print_err("allocation impossible for pthread_t");
	    exit(EXIT_FAILURE);
	}

	w->params =  malloc(w->inproc_length * sizeof (ssm_params_worker_inproc_t));
	if(w->params == NULL){
	    ssm_print_err("allocation impossible for ssm_params_worker_inproc_t");
	    exit(EXIT_FAILURE);
	}

	int J_chunk = fitness->J / w->inproc_length;
	for(i=0; i<w->inproc_length; i++){
	    w->params[i].id = opts->id;
	    w->params[i].context = w->context;
	    w->params[i].wopts = wopts;
	    w->params[i].J_chunk = J_chunk;
	    w->params[i].data = data;
	    w->params[i].J_par = J_par;
	    w->params[i].D_J_X = D_J_X;
	    w->params[i].calc = calc[i];
	    w->params[i].nav = nav;
	    w->params[i].fitness = fitness;
	    w->params[i].f_pred = f_pred;

	    pthread_create(&(w->workers[i]), NULL, ssm_worker_inproc, (void*) &(w->params[i]));
	}

	//wait that all worker are connected
	for (i = 0; i < w->inproc_length; i++) {
	    zmq_recv(w->receiver, &id, sizeof (int), 0);
	}    
    }

    return w;
}


void ssm_workers_stop(ssm_workers_t *workers)
{
    int i;

    if(workers->flag_tcp || workers->inproc_length >1){
        zmq_send (workers->controller, "KILL", 5, 0);
        zmq_close (workers->sender);
        zmq_close (workers->receiver);
        zmq_close (workers->controller);

	if(!workers->flag_tcp){
	    for(i = 0; i < workers->inproc_length; i++){
		pthread_join(workers->workers[i], NULL);
	    }
	}

        free(workers->workers);
        free(workers->params);
        zmq_ctx_destroy (workers->context);
    }

    free(workers);
}
