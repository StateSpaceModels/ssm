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
    int j, n, np1, t0, t1;

    ssm_options_t *opts = ssm_options_new();
    ssm_load_options(opts, SSM_SMC, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(nav);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t **calc = ssm_N_calc_new(jdata, nav->states_sv->length + nav->states_inc->length + nav->states_diff->length, ssm_step_ode, NULL, nav, data, fitness, opts);
    ssm_X_t ***D_J_X = ssm_D_J_X_new(data, fitness, nav->states_sv->length + nav->states_inc->length + nav->states_diff->length, opts);
    ssm_X_t ***D_J_X_tmp = ssm_D_J_X_new(data, fitness, nav->states_sv->length + nav->states_inc->length + nav->states_diff->length, opts);

    json_decref(jdata);

    ssm_input2par(par, input, calc[0], nav);
    ssm_par2X(D_J_X[0][0], par, calc[0], nav);
    for(j=1; j<fitness->J; j++){
        ssm_X_copy(D_J_X[0][j], D_J_X[0][0]);
    }

    for(j=0; j<fitness->J; j++) {
        fitness->cum_status[j] = SSM_SUCCESS;
    }

    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    for(n=0; n<data->n_obs; n++) {
        np1 = n+1;
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;

        //we are going to overwrite the content of the [np1] pointer: initialise it with values from [n]
        for(j=0;j<fitness->J;j++) {
            ssm_X_copy(D_J_X[np1][j], D_J_X[n][j]);
        }

        for(j=0;j<fitness->J;j++) {
            ssm_X_reset_inc(D_J_X[np1][j], data->rows[n]);
            fitness->cum_status[j] |= (*f_pred)(D_J_X[np1][j], t0, t1, par, nav, calc[0]);

            if(data->rows[n]->ts_nonan_length) {
                fitness->weights[j] = (fitness->cum_status[j] == SSM_SUCCESS) ?  exp(ssm_log_likelihood(data->rows[n], t1, D_J_X[np1][j], par, calc[0], nav)) : 0.0;
                fitness->cum_status[j] = SSM_SUCCESS;
            }
        }

        if(data->rows[n]->ts_nonan_length) {
            if(ssm_weight(fitness, data->rows[n], n)) {
                ssm_systematic_sampling(fitness, calc[0], n);
            }
            ssm_resample_X(fitness, &(D_J_X[np1]), &(D_J_X_tmp[np1]), n);
        }
    }

    json_decref(jparameters);

    ssm_D_J_X_free(D_J_X, data, fitness);
    ssm_N_calc_free(calc, nav);

    ssm_data_free(data);
    ssm_nav_free(nav);

    ssm_input_free(input);
    ssm_par_free(par);
   
    ssm_fitness_free(fitness);

    return 0;
}
