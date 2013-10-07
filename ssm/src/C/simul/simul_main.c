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
    int i, j, n, t0, t1;

    ssm_options_t *opts = ssm_options_new();
    ssm_load_options(opts, SSM_SIMUL, argc, argv);
    
    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    json_t *jresource = json_object_get(jparameters, "jresource");
    json_t *jprediction = NULL;
    json_t *jvalues = NULL;
    json_t *jforced = NULL;
    for(i=0; i< json_array_size(jresource); i++){
        json_t *el = json_array_get(jresource, i);
        const char* name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "prediction") == 0) {
            jprediction = json_object_get(el, "data");
            break;
        }
    }

    if(jprediction){
	jforced = json_object_get(jprediction, "forced");
	jvalues = json_object_get(jprediction, "values");
	opts->J = json_array_size(jvalues);
    }

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav, data, fitness, opts);
    ssm_X_t **J_X = ssm_J_X_new(fitness, nav, opts);
    ssm_hat_t *hat = ssm_hat_new(nav);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);

    ssm_par_t **J_par = malloc(fitness->J * sizeof (ssm_par_t *));
    if(J_par == NULL) {
        ssm_print_err("Allocation impossible for ssm_par_t *");
        exit(EXIT_FAILURE);
    }
    for(j=0; j<fitness->J; j++) {
	J_par[j] = ssm_par_new(input, calc[0], nav);
	//force parameters from values
	json_t *jval = json_array_get(jvalues, j);
	if(key in forced){
	    gsl_vector_set(J_par[j], ,);
	} else {
	    it->p[i]->f_user2par(gsl_vector_get(input, it->p[i]->offset), input, calc)
	}

	ssm_par2X(J_X[j], J_par[j], calc[0], nav);
    }

    
    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);
    
    for(j=0; j<fitness->J; j++) {
        fitness->cum_status[j] = SSM_SUCCESS;
    }

    for(n=0; n<data->length; n++) {
	t0 = (n) ? data->rows[n-1]->time: 0;
	t1 = data->rows[n]->time;

        for(j=0;j<fitness->J;j++) {
	    ssm_X_reset_inc(J_X[j], data->rows[n], nav);
	    fitness->cum_status[j] |= f_pred(J_X[j], t0, t1, J_par[j], nav, calc[0]);
	}

	if (nav->print & SSM_PRINT_HAT) {
	    ssm_hat_eval(hat, J_X, J_par, nav, calc[0], NULL, t1, 0);
	    ssm_print_hat(stdout, hat, nav, data->rows[n]);
        }

	if (nav->print & SSM_PRINT_X) {
	    for(j=0; j<fitness->J; j++) {
		ssm_print_X(stdout, J_X[j], J_par, nav, calc[0], data->rows[n], j);
	    }
	}
    }
    
    json_decref(jparameters);

    for(j=0; j<fitness->J; j++) {
        ssm_par_free(J_par[j]);
    }
    free(J_par);

    ssm_J_X_free(X, fitness);
    ssm_N_calc_free(calc, nav);
    ssm_data_free(data);
    ssm_nav_free(nav);
    ssm_fitness_free(fitness);
    ssm_hat_free(hat);

    ssm_input_free(input);

    return 0;
}
