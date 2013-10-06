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

#define LARGEST_SUM_OF_SQUARE 1e20
#define SMALLEST_LOG_LIKE -1e20 

struct s_simplex
{
    ssm_data_t *data;
    ssm_nav_t *nav;
    ssm_calc_t *calc;
    ssm_input_t *input; 
    ssm_par_t *par; 
    ssm_X_t *X;
    ssm_fitness_t *fitness;

    int flag_prior;
    int flag_least_squares;
};


/**
 * function to **minimize**
 */
static double f_simplex(const gsl_vector *theta, void *params)
{
    unsigned int n, t0, t1;
    double fitness;

    struct s_simplex *p = (struct s_simplex *) params;
    ssm_data_t *data = p->data;
    ssm_nav_t *nav = p->nav;
    ssm_calc_t *calc = p->calc;
    ssm_input_t *input = p->input; 
    ssm_par_t *par = p->par; 
    ssm_X_t *X = p->X;
    ssm_fitness_t *fit = p->fitness;
    int flag_prior = p->flag_prior;
    int flag_least_squares = p->flag_least_squares;

    ssm_err_code_t status = SSM_SUCCESS;

    if(ssm_check_ic(par, calc) != SSM_SUCCESS){

	if (!(nav->print & SSM_QUIET)) {
	    ssm_print_warning("constraints on initial conditions have not been respected: assigning bad fitness");
	}
	fitness = (flag_least_squares) ? LARGEST_SUM_OF_SQUARE: SMALLEST_LOG_LIKE;

    } else {

	ssm_theta2input(input, (gsl_vector *) theta, nav);
	ssm_input2par(par, input, calc, nav);
	ssm_par2X(X, par, calc, nav);
	X->dt = X->dt0;

	fitness=0.0;

	for(n=0; n<data->n_obs; n++) {
	    t0 = (n) ? data->rows[n-1]->time: 0;
	    t1 = data->rows[n]->time;

            ssm_X_reset_inc(X, data->rows[n], nav);
	    status |= ssm_f_prediction_ode(X, t0, t1, par, nav, calc);

	    if(data->rows[n]->ts_nonan_length) {
		if (flag_least_squares) {
		    fitness += ssm_log_likelihood(data->rows[n], X, par, calc, nav, fit);
		} else {
		    fitness += ssm_sum_square(data->rows[n], X, par, calc, nav, fit);
		}
	    }
        }

	if( status != SSM_SUCCESS ){
	    fitness = (flag_least_squares) ? LARGEST_SUM_OF_SQUARE: SMALLEST_LOG_LIKE;
	    if (!(nav->print & SSM_QUIET)) {
		ssm_print_warning("warning: something went wrong");
	    }
	}
	
	if (flag_prior && !flag_least_squares) {
	    double log_prob_prior_value;
	    ssm_err_code_t rc = ssm_log_prob_prior(&log_prob_prior_value, (gsl_vector *) theta, nav, fit);
	    if(rc != SSM_SUCCESS && !(nav->print & SSM_QUIET)){
		ssm_print_warning("error log_prob_prior computation");
	    }
	    fitness += log_prob_prior_value;
	}

    }

    if(!flag_least_squares) {
        fitness = -fitness; //GSL simplex algo minimizes so we multiply by -1 in case of log likelihood
    }

    return fitness;    
}


int main(int argc, char *argv[])
{
    ssm_options_t *opts = ssm_options_new();
    ssm_load_options(opts, SSM_SIMPLEX, argc, argv);

    json_t *jparameters = ssm_load_json_stream(stdin);
    json_t *jdata = ssm_load_data(opts);

    ssm_nav_t *nav = ssm_nav_new(jparameters, opts);
    ssm_data_t *data = ssm_data_new(jdata, nav, opts);
    ssm_fitness_t *fitness = ssm_fitness_new(data, opts);
    ssm_calc_t *calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);
    ssm_X_t *X = ssm_X_new(nav, opts);

    json_decref(jdata);

    ssm_input_t *input = ssm_input_new(jparameters, nav);
    ssm_par_t *par = ssm_par_new(input, calc, nav);
    ssm_theta_t *theta = ssm_theta_new(input, nav);   
    ssm_var_t *var = ssm_var_new(jparameters, nav);

    int n_iter = opts->n_iter;
    double size_stop = opts->size_stop;
    
    struct s_simplex params = {data, nav, calc, input, par, X, fitness, opts->flag_prior, opts->flag_least_squares};

    if (n_iter == 0 && (nav->print & SSM_PRINT_TRACE)) {
        //simply return the sum of square or the log likelihood (can be used to do slices especially with least square where smc can't be used'...)
	ssm_print_trace(stdout, theta, nav, f_simplex(theta, &params), 0);
    } else {        
	ssm_simplex(theta, var, &params, &f_simplex, nav, size_stop, n_iter);
    }

    json_decref(jparameters);

    ssm_X_free(X);
    ssm_calc_free(calc, nav);
    ssm_data_free(data);
    ssm_nav_free(nav);
    ssm_fitness_free(fitness);

    ssm_input_free(input);
    ssm_par_free(par);
    ssm_theta_free(theta);
    ssm_var_free(var);

    return 0;
}
