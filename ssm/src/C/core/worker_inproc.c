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
    int compute_fitness = p->compute_fitness;
    int is_J_par = p->is_J_par;
    int is_D_J_X = p->is_D_J_X;
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
    int *j_par = (is_J_par) ? &j: &_zero;
    int *n_X = (is_D_J_X) ? &n: &_zero;

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

            zmq_recv(receiver, &the_id, sizeof (int), 0);
            zmq_recv(receiver, &n, sizeof (int), 0);

            t0 = (n) ? data->rows[n-1]->time: 0;
            t1 = data->rows[n]->time;

            int J_start = the_id * J_chunk;
            int J_end = (the_id+1 == calc->threads_length) ? fitness->J : (the_id+1)*J_chunk;

            for(j=J_start; j<J_end; j++ ){
                ssm_X_reset_inc(D_J_X[*n_X][j], data->rows[n], nav);
                fitness->cum_status[j] |= (*f_pred)(D_J_X[*n_X][j], t0, t1, J_par[*j_par], nav, calc);

                if(compute_fitness && data->rows[n]->ts_nonan_length) {
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
