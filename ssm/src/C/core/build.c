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

void ssm_input_free(ssm_input_t *input)
{
    gsl_vector_free(input);
}

ssm_par_t *ssm_par_new(ssm_nav_t *nav)
{
    return gsl_vector_calloc(nav->par_all->length);
}

void ssm_par_free(ssm_par_t *par)
{
    gsl_vector_free(par);
}

ssm_theta_t *ssm_theta_new(ssm_nav_t *nav)
{
    return gsl_vector_calloc(nav->theta_all->length);
}

void ssm_theta_free(ssm_theta_t *theta)
{
    gsl_vector_free(theta);
}

ssm_var_t *ssm_var_new(ssm_nav_t *nav, json_t *jparameters)
{
    gsl_matrix *m = gsl_matrix_calloc(nav->theta_all->length, nav->theta_all->length);


    int i,j, index;
    ssm_it_parameters_t *it = nav->theta_all;
    json_t *resource = json_object_get(jparameters, "resource");

    for(index=0; index< json_array_size(resource); index++){
        json_t *el = json_array_get(resource, index);

        const char* name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "covariance") == 0) {

            json_t *values = json_object_get(el, "data");

            for(i=0; i<it->length; i++){
                for(j=0; j<it->length; j++){
                    json_t *jcov_i = json_object_get(values, it->p[i]->name);
                    if(jcov_i){
                        json_t *jcov_ij = json_object_get(jcov_i, it->p[j]->name);
                        if(jcov_ij){
                            if(json_is_number(jcov_ij)) {
                                gsl_matrix_set(m, i, j, json_number_value(jcov_ij));
                            } else {
                                char str[SSM_STR_BUFFSIZE];
                                sprintf(str, "error: parameters.covariance.%s.%s is not a number\n", it->p[i]->name, it->p[j]->name);
                                ssm_print_err(str);
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            }

            break;
        }
    }

    return m;
}

void ssm_var_free(ssm_var_t *var)
{
    gsl_matrix_free(var);
}


ssm_it_states_t *_ssm_it_states_new(int length)
{
    ssm_it_states_t *it = malloc(sizeof (ssm_it_states_t));
    if (it == NULL) {
        ssm_print_err("Allocation impossible for ssm_it_states_t *");
        exit(EXIT_FAILURE);
    }


    it->length = length;
    if(length){
        it->p = malloc(length * sizeof (ssm_state_t *));
        if (it->p == NULL) {
            ssm_print_err("Allocation impossible for ssm_it_states_t *");
            exit(EXIT_FAILURE);
        }
    }
    return it;
}

void _ssm_it_states_free(ssm_it_states_t *it)
{
    free(it->p);
    free(it);
}


ssm_it_parameters_t *_ssm_it_parameters_new(int length)
{
    ssm_it_parameters_t *it = malloc(sizeof (ssm_it_parameters_t));
    if (it == NULL) {
        ssm_print_err("Allocation impossible for ssm_it_parameters_t *");
        exit(EXIT_FAILURE);
    }

    it->length = length;
    if(length){
        it->p = malloc(length * sizeof (ssm_parameter_t *));
        if (it->p == NULL) {
            ssm_print_err("Allocation impossible for ssm_it_parameters_t *");
            exit(EXIT_FAILURE);
        }
    }
    return it;
}


void _ssm_it_parameters_free(ssm_it_parameters_t *it)
{
    if(it->length){
        free(it->p);
    }
    free(it);
}


ssm_nav_t *ssm_nav_new(json_t *jparameters, ssm_options_t *opts)
{
    ssm_nav_t *nav = malloc(sizeof (ssm_nav_t));
    if (nav == NULL) {
        ssm_print_err("Allocation impossible for ssm_nav_t *");
        exit(EXIT_FAILURE);
    }

    nav->implementation = opts->implementation;
    nav->noises_off = opts->noises_off;
    nav->print = opts->print;

    nav->parameters = ssm_parameters_new(&nav->parameters_length);
    nav->states = ssm_states_new(&nav->states_length, nav->parameters);
    nav->observed = ssm_observed_new(&nav->observed_length);

    nav->states_sv = ssm_it_states_sv_new(nav->states);
    nav->states_remainders = ssm_it_states_remainders_new(nav->states);
    nav->states_inc = ssm_it_states_inc_new(nav->states);
    nav->states_diff = ssm_it_states_diff_new(nav->states);

    nav->par_all = ssm_it_parameters_all_new(nav->parameters);
    nav->par_noise = ssm_it_parameters_noise_new(nav->parameters);
    nav->par_vol = ssm_it_parameters_vol_new(nav->parameters);
    nav->par_icsv = ssm_it_parameters_icsv_new(nav->parameters);
    nav->par_icdiff = ssm_it_parameters_icdiff_new(nav->parameters);

    //theta: we over-allocate the iterators
    nav->theta_all = _ssm_it_parameters_new(nav->par_all->length);
    nav->theta_no_icsv_no_icdiff = _ssm_it_parameters_new(nav->par_all->length);
    nav->theta_icsv_icdiff = _ssm_it_parameters_new(nav->par_all->length);

    //json_t jparameters with diagonal covariance term to 0.0 won't be infered: re-compute length and content
    nav->theta_all->length = 0;
    nav->theta_no_icsv_no_icdiff = 0;
    nav->theta_icsv_icdiff->length = 0;

    int index, i;
    json_t *resource = json_object_get(jparameters, "resource");

    for(index=0; index< json_array_size(resource); index++){
        json_t *el = json_array_get(resource, index);

        const char* name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "covariance") == 0) {

            json_t *values = json_object_get(el, "data");

            //for all the parameters: if covariance term and covariance term >0.0, fill theta_*
            for(i=0; i<nav->par_all->length; i++){
                json_t *jcov_i = json_object_get(values, nav->par_all->p[i]->name);
                if(jcov_i){
                    json_t *jcov_ij = json_object_get(jcov_i, nav->par_all->p[i]->name);
                    if(jcov_ij){
                        if(!json_is_number(jcov_ij)) {
                            char str[SSM_STR_BUFFSIZE];
                            sprintf(str, "error: parameters.covariance.%s.%s is not a number\n", nav->par_all->p[i]->name, nav->par_all->p[i]->name);
                            ssm_print_err(str);
                            exit(EXIT_FAILURE);
                        }

                        if(json_number_value(jcov_ij) > 0.0){

                            if( ssm_in_par(nav->par_noise, nav->par_all->p[i]->name) ) {
                                if(!(nav->noises_off & SSM_NO_WHITE_NOISE)){
                                    nav->theta_all->p[nav->theta_all->length] = nav->par_all->p[i];
                                    nav->theta_all->length += 1;
                                }
                            } else {
                                nav->theta_all->p[nav->theta_all->length] = nav->par_all->p[i];
                                nav->theta_all->length += 1;
                            }

                            int in_icsv = ssm_in_par(nav->par_icsv, nav->par_all->p[i]->name);
                            int in_icdiff = ssm_in_par(nav->par_icdiff, nav->par_all->p[i]->name);
                            if(!in_icsv && !in_icdiff){
                                nav->theta_no_icsv_no_icdiff->p[nav->theta_no_icsv_no_icdiff->length] = nav->par_all->p[i];
                                nav->theta_no_icsv_no_icdiff->length += 1;
                            } else if (in_icsv || in_icdiff){
                                if(in_icdiff){
                                    if(!(nav->noises_off & SSM_NO_DIFF)){
                                        nav->theta_icsv_icdiff->p[nav->theta_icsv_icdiff->length] = nav->par_all->p[i];
                                        nav->theta_icsv_icdiff->length += 1;
                                    }
                                } else {
                                    nav->theta_icsv_icdiff->p[nav->theta_icsv_icdiff->length] = nav->par_all->p[i];
                                    nav->theta_icsv_icdiff->length += 1;
                                }
                            }
                        }
                    }
                }
            }

            break;
        }
    }

    return nav;
}


ssm_data_t *ssm_data_new(json_t *jdata, ssm_nav_t *nav, ssm_options_t *opts)
{
    char str[SSM_STR_BUFFSIZE];
    int i, j;

    ssm_data_t *data = malloc(sizeof (ssm_data_t));
    if (data==NULL) {
        ssm_print_err("Allocation impossible for ssm_data_t");
        exit(EXIT_FAILURE);
    }

    data->dates_t0 = ssm_load_jc1_new(jdata, "starts");

    json_t *jdata_data = json_object_get(jdata, "data");

    data->length = json_array_size(jdata_data);
    data->ts_length = nav->observed_length;

    ssm_row_t **rows = malloc(data->length * sizeof (ssm_row_t *));
    if (rows==NULL) {
        ssm_print_err("Allocation impossible for ssm_data_row_t **");
        exit(EXIT_FAILURE);
    }

    data->length_nonan = 0;
    data->ind_nonan = ssm_u1_new(data->length);

    for (i=0; i< data->length; i++){
        rows[i] = malloc(sizeof (data->length * sizeof (ssm_row_t)));
        if (rows[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_data_row_t *");
            exit(EXIT_FAILURE);
        }

        json_t *jrow = json_array_get(jdata_data, i);

        json_t *jdate = json_object_get(jrow, "date");
        if(json_is_string(jdate)) {
            rows[i]->date = strdup(json_string_value(jdate));
        } else {
            sprintf(str, "error: data[%d].date is not a string\n", i);
            ssm_print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jtime = json_object_get(jrow, "time");
        if(json_is_number(jtime)) {
            rows[i]->time = (unsigned int) json_integer_value(jtime);
        } else {
            sprintf(str, "error: data[%d].time is not an integer\n", i);
            ssm_print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jobserved = json_object_get(jrow, "observed");
        rows[i]->ts_nonan_length = json_array_size(jobserved);
        rows[i]->observed = malloc(rows[i]->ts_nonan_length * sizeof (ssm_observed_t *));
        if (rows[i]->observed == NULL) {
            ssm_print_err("Allocation impossible for ssm_data_row_t.observed");
            exit(EXIT_FAILURE);
        }

        for(j=0; j<rows[i]->ts_nonan_length; j++){

            json_t *jobserved_j = json_array_get(jobserved, j);
            if(json_is_number(jobserved_j)) {
                int id = json_integer_value(jobserved_j);
                rows[i]->observed[j] = nav->observed[id];
            } else {
                sprintf(str, "error: data[%d].observed[%d] is not an integer\n", i, j);
                ssm_print_err(str);
                exit(EXIT_FAILURE);
            }
        }

        rows[i]->values = ssm_load_jd1_new(jrow, "values");

        json_t *jreset = json_object_get(jrow, "reset");
        rows[i]->states_reset_length = json_array_size(jreset);
        rows[i]->states_reset = malloc(rows[i]->states_reset_length * sizeof (ssm_state_t *));
        if (rows[i]->states_reset == NULL) {
            ssm_print_err("Allocation impossible for ssm_data_row_t.states_reset");
            exit(EXIT_FAILURE);
        }

        for(j=0; j<rows[i]->states_reset_length; j++){

            json_t *jreset_j = json_array_get(jreset, j);
            if(json_is_number(jreset_j)) {
                int id = json_integer_value(jreset_j);
                rows[i]->states_reset[j] = nav->states[id];
            } else {
                sprintf(str, "error: data[%d].reset[%d] is not an integer\n", i, j);
                ssm_print_err(str);
                exit(EXIT_FAILURE);
            }
        }

        if(rows[i]->ts_nonan_length){
            data->ind_nonan[data->length_nonan] = i;
            data->length_nonan += 1;
        }
    }

    data->rows = rows;

    //n_obs
    if(opts->n_obs >= 0){
        data->n_obs = (opts->n_obs < data->length) ? opts->n_obs : data->length;
    } else {
        data->n_obs = data->length;
    }

    return data;
}

ssm_calc_t *ssm_calc_new(json_t *jdata, int dim_ode, int (*func_step_ode) (double t, const double y[], double dydt[], void * params), int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params), ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, int thread_id, unsigned long int seed, ssm_options_t *opts)
{
    ssm_calc_t *calc = malloc(sizeof (ssm_calc_t));
    if (calc==NULL) {
        char str[SSM_STR_BUFFSIZE];
        snprintf(str, SSM_STR_BUFFSIZE, "Allocation impossible for ssm_calc_t (thread_id: %d)", thread_id);
        ssm_print_err(str);
        exit(EXIT_FAILURE);
    }

    /***********/
    /* threads */
    /***********/

    calc->threads_length = ssm_sanitize_n_threads(opts->n_thread, fitness);
    calc->thread_id = thread_id;

    /******************/
    /* random numbers */
    /******************/

    //  random number generator and parallel MC simulations:
    //
    //  idea using one different seed per thread but is it realy uncorelated ???
    //  Should I go through the trouble of changing from GSL to SPRNG????
    //  answer:
    //  I would recommend using ranlxd.  The seeds should give 2^31
    //  effectively independent streams of length 10^171.  A discussion of the
    //  seeding procedure can be found in the file notes.ps at
    //  http://www.briangough.ukfsn.org/ranlux_2.2/
    //  --
    //  Brian Gough
    //
    //  => we create as many rng as parallel threads *but* note that for
    //  the operations not prarallelized, we always use
    //  cacl[0].randgsl


    const gsl_rng_type *Type;
    if (calc->threads_length == 1){ //we don't need a rng supporting parallel computing, we use mt19937 that is way faster than ranlxs0 (1754 k ints/sec vs 565 k ints/sec)
        Type = gsl_rng_mt19937; /*MT19937 generator of Makoto Matsumoto and Takuji Nishimura*/
    } else {
        Type = gsl_rng_ranlxs0; //gsl_rng_ranlxs2 is better than gsl_rng_ranlxs0 but 2 times slower
    }

    calc->randgsl = gsl_rng_alloc(Type);
    gsl_rng_set(calc->randgsl, seed + thread_id);

    /*******************/
    /* implementations */
    /*******************/

    if (nav->implementation == SSM_ODE || nav->implementation == SSM_EKF){

        calc->T = gsl_odeiv2_step_rkf45;
        calc->control = gsl_odeiv2_control_y_new(opts->eps_abs, opts->eps_rel);
        calc->step = gsl_odeiv2_step_alloc(calc->T, dim_ode);
        calc->evolve = gsl_odeiv2_evolve_alloc(dim_ode);
        (calc->sys).function = func_step_ode;
        (calc->sys).jacobian = jacobian;
        (calc->sys).dimension= dim_ode;
        (calc->sys).params= calc;

        calc->yerr = ssm_d1_new(dim_ode);

        if(nav->implementation == SSM_EKF){
            int n_s = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
            int n_o = nav->observed_length;
            calc->_pred_error = gsl_vector_calloc(n_o);
            calc->_St = gsl_matrix_calloc(n_o, n_o);
            calc->_Stm1 = gsl_matrix_calloc(n_o, n_o);
            calc->_Rt = gsl_matrix_calloc(n_o, n_o);
            calc->_Ht = gsl_matrix_calloc(n_s, n_o);
            calc->_Kt = gsl_matrix_calloc(n_s, n_o);
            calc->_Tmp_N_SV_N_TS = gsl_matrix_calloc(n_s, n_o);
            calc->_Tmp_N_TS_N_SV = gsl_matrix_calloc(n_o, n_s);
            calc->_Q = gsl_matrix_calloc(n_s, n_s);
            calc->_FtCt = gsl_matrix_calloc(n_s, n_s);
            calc->_Ft = gsl_matrix_calloc(n_s, n_s);
        }

    } else if (nav->implementation == SSM_SDE){
        calc->y_pred = ssm_d1_new(dim_ode);
    } else if (nav->implementation == SSM_PSR){
        ssm_alloc_psr(calc);
    }

    /**************************/
    /* multi-threaded sorting */
    /**************************/

    calc->to_be_sorted = ssm_d1_new(fitness->J);
    calc->index_sorted = ssm_st1_new(fitness->J);

    /**************/
    /* covariates */
    /**************/

    json_t *jcovariates = json_object_get(jdata, "covariates");
    calc->covariates_length = json_array_size(jcovariates);

    calc->acc = malloc(calc->covariates_length * sizeof(gsl_interp_accel *));
    if (calc->acc == NULL) {
        char str[SSM_STR_BUFFSIZE];
        snprintf(str, SSM_STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        ssm_print_err(str);
        exit(EXIT_FAILURE);
    }

    calc->spline = malloc(calc->covariates_length * sizeof(gsl_spline *));
    if (calc->spline == NULL) {
        char str[SSM_STR_BUFFSIZE];
        snprintf(str, SSM_STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        ssm_print_err(str);
        exit(EXIT_FAILURE);
    }

    const gsl_interp_type *my_gsl_interp_type = ssm_str_to_interp_type(opts->interpolator);

    int k, z;

    //TODO init from opts !!!!
    double freeze_forcing = -1.0; // the time (in days) to freeze (i.e only take metadata from this time) (ignored if freeze_forcing < 0.0)
    double t_max = -1.0; //t_max the highest possible time in days when interpolated metadata will be requested (-1 to default to last point of metadata). t_

    for (k=0; k< calc->covariates_length; k++) {
        json_t *jcovariate = json_array_get(jcovariates, k);

        double *x = ssm_load_jd1_new(jcovariate, "x");
        double *y = ssm_load_jd1_new(jcovariate, "y");
        int size = json_array_size(json_object_get(jcovariate, "x"));

        if((freeze_forcing < 0.0) && (t_max > x[size-1])){ //no freeze but t_max > x[size-1] repeat last value
            int prev_size = size ;
            size += ((int) t_max - x[prev_size-1]) ;

            double *tmp_x = realloc(x, size * sizeof (double) );
            if ( tmp_x == NULL ) {
                ssm_print_err("Reallocation impossible"); free(x); exit(EXIT_FAILURE);
            } else {
                x = tmp_x;
            }

            double *tmp_y = realloc(y, size * sizeof (double) );
            if ( tmp_y == NULL ) {
                ssm_print_err("Reallocation impossible"); free(y); exit(EXIT_FAILURE);
            } else {
                y = tmp_y;
            }

            //repeat last value
            double xlast = x[prev_size-1];
            for(z = prev_size;  z < size ; z++ ){
                x[z] = xlast + z;
                y[z] = y[prev_size - 1];
            }
        }

        if( (freeze_forcing>=0.0) || (size == 1) ){ //only 1 value: make it 2
            double x_all[2];
            x_all[0] = x[0];
            x_all[1] = GSL_MAX(GSL_MAX(t_max, ((data->n_obs>=1) ? (double) data->rows[data->n_obs-1]->time: 0.0)), x[size-1]);

            double y_all[2];
            y_all[0] = (size == 1) ? y[0]: gsl_spline_eval(calc->spline[k], GSL_MIN(freeze_forcing, x[size-1]), calc->acc[k]); //interpolate y for time freeze_forcing requested (if possible)
            y_all[1] = y_all[0];

            calc->acc[k] = gsl_interp_accel_alloc ();
            calc->spline[k]  = gsl_spline_alloc (gsl_interp_linear, 2);
            gsl_spline_init (calc->spline[k], x_all, y_all, 2);

        } else if (size >= gsl_interp_type_min_size(my_gsl_interp_type)) {

            calc->acc[k] = gsl_interp_accel_alloc ();
            calc->spline[k]  = gsl_spline_alloc(my_gsl_interp_type, size);
            gsl_spline_init (calc->spline[k], x, y, size);

        } else {

            ssm_print_warning("insufficient data points for required metadata interpolator, switching to linear");
            calc->acc[k] = gsl_interp_accel_alloc ();
            calc->spline[k] = gsl_spline_alloc (gsl_interp_linear, size);
            gsl_spline_init(calc->spline[k], x, y, size);

        }

        free(x);
        free(y);
    }

    return calc;
}


ssm_calc_t **ssm_N_calc_new(json_t *jdata, int dim_ode, int (*func_step_ode) (double t, const double y[], double dydt[], void * params), int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params), ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, ssm_options_t *opts)
{
    int i;
    int n_threads = ssm_sanitize_n_threads(opts->n_thread, fitness);

    ssm_calc_t **calc = malloc(n_threads * sizeof (ssm_calc_t *));
    if (calc==NULL) {
        ssm_print_err("Allocation impossible for ssm_calc_t **");
        exit(EXIT_FAILURE);
    }

    unsigned long int seed;
    if(opts->flag_seed_time){
        seed = (unsigned) time(NULL);
    } else{
        seed=2;
    }
    seed += opts->id; /*we ensure uniqueness of seed in case of parrallel runs*/

    for (i=0; i< n_threads; i++) {
        calc[i] = ssm_calc_new(jdata, dim_ode, func_step_ode, jacobian, nav, data, fitness, i, seed, opts);
    }

    return calc;
}
