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
    int n, t0, t1;
    ssm_options_t *opts = ssm_options_new();
    ssm_options_load(opts, SSM_KALMAN, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t *calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);
    ssm_X_t *X = ssm_X_new(nav, opts);
    ssm_hat_t *hat = ssm_hat_new(nav);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc, nav);
    ssm_theta_t *theta = ssm_theta_new(input, nav);

    int flag_prior = opts->flag_prior;
    int flag_no_filter = opts->flag_no_filter;

    ssm_f_pred_t f_pred = ssm_get_f_pred(nav);

    fitness->cum_status[0] = SSM_SUCCESS;

    ssm_par2X(X, par, calc, nav);

    for(n=0; n<data->n_obs; n++) {
        t0 = (n) ? data->rows[n-1]->time: 0;
        t1 = data->rows[n]->time;

        ssm_X_reset_inc(X, data->rows[n], nav);
        fitness->cum_status[0] |= (*f_pred)(X, t0, t1, par, nav, calc);

        if (nav->print & SSM_PRINT_DIAG) {
            ssm_print_pred_res(nav->diag, &X, par, nav, calc, data, data->rows[n], fitness);
        }
	
        if(!flag_no_filter && data->rows[n]->ts_nonan_length && (fitness->cum_status[0] == SSM_SUCCESS)) {
            fitness->cum_status[0] |= ssm_kalman_update(fitness, X, data->rows[n], t1, par, calc, nav);
	    if( fitness->cum_status[0] & SSM_ERR_REM_SV ) { 
		fitness->log_like =  GSL_NEGINF;
		if(nav->print & SSM_PRINT_WARNING) {
		    ssm_print_warning("error: negative state variable or remainder");
		}
	    }
        }

        if (nav->print & SSM_PRINT_HAT) {
            ssm_hat_eval(hat, &X, &par, nav, calc, fitness, t1, 0);
            ssm_print_hat(nav->hat, hat, nav, data->rows[n]);
        }

        if (nav->print & SSM_PRINT_X) {
            ssm_print_X(nav->X, X, par, nav, calc, data->rows[n], 0);
        }
    }

    if (flag_prior) {
        double log_prob_prior_value;
        ssm_err_code_t rc = ssm_log_prob_prior(&log_prob_prior_value, theta, nav, fitness);
        if(rc != SSM_SUCCESS && (nav->print & SSM_PRINT_WARNING)){
            ssm_print_warning("error log_prob_prior computation");
        }
	fitness->summary_log_ltp = fitness->log_like + log_prob_prior_value;
    } else {
	ssm_aic(fitness, nav, fitness->log_like);
    }

    if (nav->print & SSM_PRINT_TRACE) {
        ssm_print_trace(nav->trace, theta, nav, (flag_prior)? fitness->summary_log_ltp : fitness->log_like, 0);
    }

    if (!(nav->print & SSM_PRINT_LOG)) {
	ssm_pipe_theta(stdout, jparameters, theta, NULL, fitness, nav, opts);
    } else {
	char str[SSM_STR_BUFFSIZE];
	snprintf(str, SSM_STR_BUFFSIZE, "log(like%s): %g", (flag_prior)? "*prior": "", fitness->log_like);
	ssm_print_log(str);
    }


    json_decref(jparameters);

    ssm_X_free(X);
    ssm_calc_free(calc, nav);
    ssm_hat_free(hat);

    ssm_data_free(data);
    ssm_nav_free(nav);

    ssm_input_free(input);
    ssm_par_free(par);
    ssm_theta_free(theta);

    ssm_fitness_free(fitness);

    return 0;
}
