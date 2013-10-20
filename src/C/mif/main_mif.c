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
    int j, n, np1, t0, t1, the_j;

    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_MIF, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav, data, fitness, opts);
    ssm_X_t **J_X = ssm_J_X_new(fitness, nav, opts);
    ssm_X_t **J_X_tmp = ssm_J_X_new(fitness, nav, opts);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_theta_t *mle = ssm_theta_new(input, nav);
    ssm_var_t *var = ssm_var_new(jparameters, nav);
    ssm_mif_scale_var(var, data, nav);

    ssm_par_t **J_par = malloc(fitness->J * sizeof (ssm_par_t *));
    if(J_par == NULL) {
        ssm_print_err("Allocation impossible for ssm_par_t *");
        exit(EXIT_FAILURE);
    }
    ssm_theta_t **J_theta = malloc(fitness->J * sizeof (ssm_theta_t *));
    if(J_theta == NULL) {
        ssm_print_err("Allocation impossible for ssm_theta_t *");
        exit(EXIT_FAILURE);
    }
    ssm_theta_t **J_theta_tmp = malloc(fitness->J * sizeof (ssm_theta_t *));
    if(J_theta_tmp == NULL) {
        ssm_print_err("Allocation impossible for ssm_theta_t *");
        exit(EXIT_FAILURE);
    }

    for(j=0; j<fitness->J; j++) {
        J_par[j] = ssm_par_new(input, calc[0], nav);
        J_theta[j] = ssm_theta_new(NULL, nav);
        J_theta_tmp[j] = ssm_theta_new(NULL, nav);
    }

    double **D_theta_bart = ssm_d2_new(data->length+1, nav->theta_all->length); //mean of theta at each time step, +1 because we keep values for every data point + initial condition
    double **D_theta_Vt = ssm_d2_new(data->length+1, nav->theta_all->length); //variance of theta at each time step

    int m, i, id;
    int n_iter = opts->n_iter;
    int flag_prior = opts->flag_prior;
    int L = (int) floor(opts->L*data->length);
    double delta;
    double cooling;

    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    ssm_workers_t *workers = ssm_workers_start(&J_X, J_par, data, calc, fitness, f_pred, nav, opts, SSM_WORKER_J_PAR | SSM_WORKER_FITNESS);

    for(m=1; m<n_iter; m++){

        fitness->log_like = 0.0;
        fitness->n_all_fail = 0;
        delta = 0;
        cooling = ssm_mif_cooling(opts, m);

        for(i=0; i<nav->theta_all->length; i++){
            D_theta_bart[0][i] = gsl_vector_get(mle, i);
            D_theta_Vt[0][i] = ( pow(opts->b * cooling, 2) * gsl_matrix_get(var, i, i) );
        }

        for(j=0; j<fitness->J; j++) {
            do{
                ssm_theta_ran(J_theta[j], mle, var, opts->b * cooling, calc[0], nav, 0);
                ssm_theta2input(input, J_theta[j], nav);
                ssm_input2par(J_par[j], input, calc[0], nav);
            } while(ssm_check_ic(J_par[j], calc[0]) != SSM_SUCCESS);

            ssm_par2X(J_X[j], J_par[j], calc[0], nav);
            J_X[j]->dt = J_X[0]->dt0;
            fitness->cum_status[j] = SSM_SUCCESS;
        }

        for(n=0; n<data->n_obs; n++) {
            np1 = n+1;
            t0 = (n) ? data->rows[n-1]->time: 0;
            t1 = data->rows[n]->time;
            delta += (t1-t0); //cumulate t1-t0 in between 2 data step where data->rows[n]->ts_nonan_length > 0


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
		    fitness->cum_status[j] |= (*f_pred)(J_X[j], t0, t1, J_par[j], nav, calc[0]);

		    if(data->rows[n]->ts_nonan_length) {
			fitness->weights[j] = (fitness->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], J_X[j], J_par[j], calc[0], nav, fitness)) : 0.0;
			fitness->cum_status[j] = SSM_SUCCESS;
		    }
		}
	    }


            if(data->rows[n]->ts_nonan_length) {
                if (flag_prior) {
                    ssm_mif_patch_like_prior(fitness->weights, fitness, J_theta, data, nav, n, L);
                }

                int some_particle_succeeded = ssm_weight(fitness, data->rows[n], nav, n);

                ssm_mif_mean_var_theta_theoretical(D_theta_bart[np1], D_theta_Vt[np1], J_theta, var, fitness, nav, delta*pow(cooling, 2));
                if (nav->print & SSM_PRINT_DIAG) {
                    ssm_mif_print_mean_var_theoretical_ess(nav->diag, D_theta_bart[np1], D_theta_Vt[np1], fitness, nav , data->rows[n], m);
                }

                if(some_particle_succeeded) {
                    ssm_systematic_sampling(fitness, calc[0], n);
                }

                ssm_mif_resample_and_mutate_theta(fitness, J_theta, J_theta_tmp, var, calc, nav, cooling*sqrt(delta), n);
                ssm_resample_X(fitness, &J_X, &J_X_tmp, n);

                delta = 0.0;
            }

            if(n == L){
                ssm_mif_fixed_lag_smoothing(mle, J_theta, fitness, nav);
            }
        }

        (m<=opts->m_switch) ? ssm_mif_update_average(mle, D_theta_bart, data, nav): ssm_mif_update_ionides(mle, var, D_theta_bart, D_theta_Vt, data, nav, opts, cooling);

        if(nav->print & SSM_PRINT_TRACE){
            ssm_print_trace(nav->trace, mle, nav, fitness->log_like, m);
        }

    }

    if (!(nav->print & SSM_PRINT_LOG)) {
	ssm_pipe_theta(stdout, jparameters, mle, NULL, nav, opts);
    }

    json_decref(jparameters);

    ssm_workers_stop(workers);

    ssm_J_X_free(J_X, fitness);
    ssm_J_X_free(J_X_tmp, fitness);
    ssm_N_calc_free(calc, nav);

    ssm_d2_free(D_theta_bart, data->length+1);
    ssm_d2_free(D_theta_Vt, data->length+1);

    ssm_data_free(data);
    ssm_nav_free(nav);
    ssm_fitness_free(fitness);

    ssm_var_free(var);
    ssm_input_free(input);
    ssm_theta_free(mle);

    for(j=0; j<fitness->J; j++) {
        ssm_par_free(J_par[j]);
        ssm_theta_free(J_theta[j]);
        ssm_theta_free(J_theta_tmp[j]);
    }
    free(J_par);
    free(J_theta);
    free(J_theta_tmp);

    return 0;
}
