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

    int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;

    ssm_X_t **D_X = ssm_D_J_X_new(data->n_data +1, 1, m + pow(m,2), dt);

    ssm_fitness_t *like = ssm_fitness_new(1, data, like_min);

    ssm_calc_t **calc = ssm_N_calc_new(n_threads, ...);

    json_decref(settings);

    ssm_par2X(D_J_X[0][0], par, calc[0], nav);

    int j, n, np1, t0,t1;

    like->cum_status[0] = SSM_SUCCESS;
    
    for(n=0; n<data->n_obs; n++) {
	np1 = n+1;
	t0 = data->times[n];
	t1 = data->times[np1];
		    
	//we are going to overwrite the content of the [np1] pointer: initialise it with values from [n]
	ssm_X_copy(D_J_X[np1][0], D_J_X[n][0]);
	ssm_X_reset_inc(D_J_X[np1][0], data->rows[n]);
	
	// Predict
	like->cum_status[0] |= (*f_pred)(D_J_X[np1][0], t0, t1, par, nav, calc[0]);

	// Update
	like->cum_status[0] |= ssm_kalman_update(D_J_X[np1][0], data->rows[n], t1, par, calc[0], nav, like)
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
