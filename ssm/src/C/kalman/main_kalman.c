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
    int j, n, np1, t0,t1;

    ssm_options_t *opts = ssm_options_new();
    ssm_load_options(opts, SSM_KALMAN, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(nav);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t *calc = ssm_calc_new(jdata, nav->states_sv_inc->length + nav->states_diff->length, ssm_step_ode, NULL, nav, data, fitness, opts);
    ssm_X_t *X = ssm_X_new(data->n_data +1, nav);
 
    json_decref(jdata);

    ssm_input2par(par, input, calc, nav);
    ssm_par2X(X, par, calc, nav);

    fitness->cum_status[0] = SSM_SUCCESS;
    
    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    for(n=0; n<data->n_obs; n++) {
	np1 = n+1;
	t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;
		    
	// Reset incidence
	ssm_X_reset_inc_and_cov(X, data->rows[n], nav);
	
	// Predict
	fitness->cum_status[0] |= (*f_pred)(X, t0, t1, par, nav, calc);

	// Update
	fitness->cum_status[0] |= ssm_kalman_update(X, data->rows[n], t1, par, calc, nav, fitness);
    }

    json_decref(parameters);

    ssm_X_free(X);
    ssm_calc_free(calc);

    ssm_data_free(data);
    ssm_nav_free(nav);
    
    ssm_input_free(input);
    ssm_par_free(par);
    
    ssm_fitness_free(fitness);

    return 0;
}
