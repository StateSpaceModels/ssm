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

int main(int argc, char *argv[]){

    json_t parameters = load_json();
    json_t *settings = load_settings(SSM_PATH_SETTINGS);

    ssm_nav_t *nav = ssm_nav_new(settings, parameters);
    ssm_data_t *data = ssm_data_new(settings, n_obs);

    
    ssm_input_t *input = ssm_input_new(parameters, nav);
    ssm_par_t *par = ssm_par_new(nav);
    ssm_input2par(par, input, calc, nav);

    ssm_X_t ***D_J_X = ssm_D_J_X_new(data->n_data +1, J, nav->states_sv->length + nav->states_inc->length + nav->states_diff->length, dt);

    ssm_fitness_t *like = ssm_fitness_new(J, data, like_min);

    ssm_calc_t **calc = ssm_N_calc_new(n_threads, ...);

    json_decref(settings);

    ssm_par2X(D_J_X[0][0], par, calc[0], nav);
    for(j=1; j<J; j++){
	ssm_X_copy(D_J_p_X[0][j], D_J_p_X[0][0]);
    }

    int j, n, np1, t0,t1;

    for(j=0; j<like->J; j++) {
	like->cum_status[j] = SSM_SUCCESS;
    }

    for(n=0; n<data->n_obs; n++) {
	np1 = n+1;
	t0 = (n) ? data->rows[n-1]->time: 0.0;
	t1 = data->rows[n]->time;

	//we are going to overwrite the content of the [np1] pointer: initialise it with values from [n]
	for(j=0;j<like->J;j++) {
	    ssm_X_copy(D_J_X[np1][j], D_J_X[n][j]);
	}

	for(j=0;j<like->J;j++) {
	    ssm_X_reset_inc(D_J_X[np1][j], data->rows[n]);
	    like->cum_status[j] |= (*f_pred)(D_J_X[np1][j], t0, t1, par, nav, calc[0]);

	    if(data->rows[n]->ts_nonan_length) {
		like->weights[j] = (like->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], t1, D_J_X[np1][j], par, calc[0], nav)) : 0.0;
		like->cum_status[j] = SSM_SUCCESS;
	    }
	}

        if(data->rows[n]->ts_nonan_length) {
	    if(ssm_weight(like, n)) {
		ssm_systematic_sampling(like, calc[0], n);
	    }
	    ssm_resample_X(like, &(D_J_X[np1]), &(D_J_X_tmp[np1]), n);
	}
    }

    json_decref(parameters);

    ssm_data_free(data);
    ssm_nav_free(nav);
    
    ssm_input_free(input);
    ssm_par_free(par);
    
    ssm_D_J_X_free(D_J_X);

    ssm_N_calc_free(calc);

    ssm_fitness_free(like);

    return 0;
}
