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

/**
 * !!MUST be construced from ssm_input_t
 */
ssm_par_t *ssm_par_new(ssm_input_t *input, ssm_calc_t *calc, ssm_nav_t *nav)
{
    ssm_it_parameters_t *it = nav->par_all;
    ssm_par_t *par = gsl_vector_calloc(it->length);
    ssm_input2par(par, input, calc, nav);

    return par;
}

void ssm_par_free(ssm_par_t *par)
{
    gsl_vector_free(par);
}

/**
 * if input is NULL return empty theta
 */
ssm_theta_t *ssm_theta_new(ssm_input_t* input, ssm_nav_t *nav)
{
    ssm_theta_t *theta = gsl_vector_calloc(nav->theta_all->length);

    if(input) {
        int i;
        ssm_parameter_t *p;

        for(i=0; i< nav->theta_all->length; i++){
            p = nav->theta_all->p[i];
            gsl_vector_set(theta, i, p->f(gsl_vector_get(input, p->offset)));
        }
    }

    return theta;
}

void ssm_theta_free(ssm_theta_t *theta)
{
    gsl_vector_free(theta);
}

ssm_var_t *ssm_var_new(json_t *jparameters, ssm_nav_t *nav)
{
    gsl_matrix *m = gsl_matrix_calloc(nav->theta_all->length, nav->theta_all->length);


    int i,j, index;
    ssm_it_parameters_t *it = nav->theta_all;
    json_t *jresource = json_object_get(jparameters, "resources");

    for(index=0; index< json_array_size(jresource); index++){
        json_t *el = json_array_get(jresource, index);

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
                                snprintf(str, SSM_STR_BUFFSIZE, "error: parameters.covariance.%s.%s is not a number\n", it->p[i]->name, it->p[j]->name);
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
    it->p = malloc(length * sizeof (ssm_state_t *));
    if (length && (it->p == NULL)) {
        ssm_print_err("Allocation impossible for ssm_it_states_t *");
        exit(EXIT_FAILURE);
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
    it->p = malloc(length * sizeof (ssm_parameter_t *));
    if (length && (it->p == NULL)) {
        ssm_print_err("Allocation impossible for ssm_it_parameters_t *");
        exit(EXIT_FAILURE);
    }

    return it;
}


void _ssm_it_parameters_free(ssm_it_parameters_t *it)
{
    free(it->p);
    free(it);
}


ssm_nav_t *ssm_nav_new(json_t *jparameters, ssm_options_t *opts)
{
    char str[SSM_STR_BUFFSIZE];

    ssm_nav_t *nav = malloc(sizeof (ssm_nav_t));
    if (nav == NULL) {
        ssm_print_err("Allocation impossible for ssm_nav_t *");
        exit(EXIT_FAILURE);
    }

    nav->implementation = opts->implementation;
    nav->noises_off = opts->noises_off;
    nav->print = opts->print;

    nav->parameters = _ssm_parameters_new(&nav->parameters_length);
    nav->states = _ssm_states_new(&nav->states_length, nav->parameters);
    nav->observed = _ssm_observed_new(&nav->observed_length);

    nav->states_sv = ssm_it_states_sv_new(nav->states);
    nav->states_remainders = ssm_it_states_remainders_new(nav->states);
    nav->states_inc = ssm_it_states_inc_new(nav->states);
    nav->states_sv_inc = ssm_it_states_sv_inc_new(nav->states);
    nav->states_diff = ssm_it_states_diff_new(nav->states);

    nav->par_all = ssm_it_parameters_all_new(nav->parameters);
    nav->par_noise = ssm_it_parameters_noise_new(nav->parameters);
    nav->par_disp = ssm_it_parameters_disp_new(nav->parameters);
    nav->par_icsv = ssm_it_parameters_icsv_new(nav->parameters);
    nav->par_icdiff = ssm_it_parameters_icdiff_new(nav->parameters);

    //theta: we over-allocate the iterators
    nav->theta_all = _ssm_it_parameters_new(nav->par_all->length);
    nav->theta_no_icsv_no_icdiff = _ssm_it_parameters_new(nav->par_all->length);
    nav->theta_icsv_icdiff = _ssm_it_parameters_new(nav->par_all->length);

    //json_t jparameters with diagonal covariance term to 0.0 won't be infered: re-compute length and content
    nav->theta_all->length = 0;
    nav->theta_no_icsv_no_icdiff->length = 0;
    nav->theta_icsv_icdiff->length = 0;


    int index, i;
    json_t *jresource = json_object_get(jparameters, "resources");

    for(index=0; index< json_array_size(jresource); index++){
        json_t *el = json_array_get(jresource, index);

        const char *name = json_string_value(json_object_get(el, "name"));
        if (strcmp(name, "covariance") == 0) {

            json_t *values = json_object_get(el, "data");

            //for all the parameters: if covariance term and covariance term >0.0, fill theta_*
            for(i=0; i<nav->par_all->length; i++){
		json_t *jcov_i = json_object_get(values, nav->par_all->p[i]->name);
                if(jcov_i){
                    json_t *jcov_ii = json_object_get(jcov_i, nav->par_all->p[i]->name);
                    if(jcov_ii){
                        if(!json_is_number(jcov_ii)) {
                            snprintf(str, SSM_STR_BUFFSIZE, "error: parameters.covariance.%s.%s is not a number\n", nav->par_all->p[i]->name, nav->par_all->p[i]->name);
                            ssm_print_err(str);
                            exit(EXIT_FAILURE);
                        }

                        if(json_number_value(jcov_ii) > 0.0){

                            if( ssm_in_par(nav->par_noise, nav->par_all->p[i]->name) ) {
                                if(!(nav->noises_off & SSM_NO_WHITE_NOISE)){
                                    nav->theta_all->p[nav->theta_all->length] = nav->par_all->p[i];
                                    nav->theta_all->p[nav->theta_all->length]->offset_theta = nav->theta_all->length;
                                    nav->theta_all->length += 1;
                                }
                            } else if( ssm_in_par(nav->par_disp, nav->par_all->p[i]->name) ) {
                                if(!(nav->noises_off & SSM_NO_DIFF)){
                                    nav->theta_all->p[nav->theta_all->length] = nav->par_all->p[i];
                                    nav->theta_all->p[nav->theta_all->length]->offset_theta = nav->theta_all->length;
                                    nav->theta_all->length += 1;
                                }
                            } else {
                                nav->theta_all->p[nav->theta_all->length] = nav->par_all->p[i];
                                nav->theta_all->p[nav->theta_all->length]->offset_theta = nav->theta_all->length;
                                nav->theta_all->length += 1;
                            }

                            int in_icsv = ssm_in_par(nav->par_icsv, nav->par_all->p[i]->name);
                            int in_icdiff =ssm_in_par(nav->par_icdiff, nav->par_all->p[i]->name);

                            if(!in_icsv && !in_icdiff){
                                nav->theta_no_icsv_no_icdiff->p[nav->theta_no_icsv_no_icdiff->length] = nav->par_all->p[i];
                                nav->theta_no_icsv_no_icdiff->length += 1;
                            } else if (in_icsv || in_icdiff){
                                if(in_icdiff){
                                    if(!(nav->noises_off & SSM_NO_DIFF)){ //diffusion is allowed
                                        nav->theta_icsv_icdiff->p[nav->theta_icsv_icdiff->length] = nav->par_all->p[i];
                                        nav->theta_icsv_icdiff->length += 1;
                                    } else {
                                        nav->theta_no_icsv_no_icdiff->p[nav->theta_no_icsv_no_icdiff->length] = nav->par_all->p[i];
                                        nav->theta_no_icsv_no_icdiff->length += 1;
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

    //files (CSV open and print headers)
    if(opts->print & SSM_PRINT_TRACE){
#if SSM_JSON
        nav->trace = stdout;
#else
        snprintf(str, SSM_STR_BUFFSIZE, "%s/trace_%d.csv", opts->root,  opts->id);
        nav->trace = fopen(str, "w");
        ssm_print_header_trace(nav->trace, nav);
#endif
    } else {
        nav->trace = NULL;
    }

    if(opts->print & SSM_PRINT_X){
#if SSM_JSON
        nav->X = stdout;
#else
snprintf(str, SSM_STR_BUFFSIZE, "%s/X_%d.csv", opts->root,  opts->id);
 nav->X = fopen(str, "w");
 ssm_print_header_X(nav->X, nav);
#endif
    } else {
        nav->X = NULL;
    }

    if(opts->print & SSM_PRINT_HAT){
#if SSM_JSON
        nav->hat = stdout;
#else
        snprintf(str, SSM_STR_BUFFSIZE, "%s/hat_%d.csv", opts->root,  opts->id);
        nav->hat = fopen(str, "w");
        ssm_print_header_hat(nav->hat, nav);
#endif
    } else {
        nav->hat = NULL;
    }

    if(opts->print & SSM_PRINT_DIAG){
#if SSM_JSON
        nav->diag = stdout;
#else
        snprintf(str, SSM_STR_BUFFSIZE, "%s/diag_%d.csv", opts->root,  opts->id);
        nav->diag = fopen(str, "w");
        if(opts->algo & (SSM_SMC | SSM_KALMAN)){
            ssm_print_header_pred_res(nav->diag, nav);
        } else if (opts->algo & (SSM_PMCMC | SSM_KMCMC)){
            ssm_print_header_ar(nav->diag);
        } else if (opts->algo & SSM_MIF){
            ssm_mif_print_header_mean_var_theoretical_ess(nav->diag, nav);
        }
#endif
    } else {
        nav->diag = NULL;
    }

    return nav;
}


void _ssm_observed_free(ssm_observed_t *observed)
{
    free(observed->name);
    free(observed);
}

void _ssm_parameter_free(ssm_parameter_t *parameter)
{
    free(parameter->name);
    free(parameter);
}

void _ssm_state_free(ssm_state_t *state)
{
    free(state->name);
    free(state);
}


void ssm_nav_free(ssm_nav_t *nav)
{
    int i;

    _ssm_it_states_free(nav->states_sv);
    _ssm_it_states_free(nav->states_remainders);
    _ssm_it_states_free(nav->states_inc);
    _ssm_it_states_free(nav->states_sv_inc);
    _ssm_it_states_free(nav->states_diff);

    _ssm_it_parameters_free(nav->par_all);
    _ssm_it_parameters_free(nav->par_noise);
    _ssm_it_parameters_free(nav->par_disp);
    _ssm_it_parameters_free(nav->par_icsv);
    _ssm_it_parameters_free(nav->par_icdiff);

    _ssm_it_parameters_free(nav->theta_all);
    _ssm_it_parameters_free(nav->theta_no_icsv_no_icdiff);
    _ssm_it_parameters_free(nav->theta_icsv_icdiff);

    for(i=0; i<nav->parameters_length; i++){
        _ssm_parameter_free(nav->parameters[i]);
    }
    free(nav->parameters);

    for(i=0; i<nav->states_length; i++){
        _ssm_state_free(nav->states[i]);
    }
    free(nav->states);

    for(i=0; i<nav->observed_length; i++){
        _ssm_observed_free(nav->observed[i]);
    }
    free(nav->observed);

#if !SSM_JSON
    if(nav->X)     fclose(nav->X);
    if(nav->hat)   fclose(nav->hat);
    if(nav->diag)  fclose(nav->diag);
    if(nav->trace) fclose(nav->trace);
#endif

    free(nav);
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

    json_t *jstart = json_object_get(jdata, "start");
    if(json_is_string(jstart)) {
        data->date_t0 = strdup(json_string_value(jstart));
    } else {
        ssm_print_err("data start is not a string");
        exit(EXIT_FAILURE);
    }


    json_t *jdata_data = json_object_get(jdata, "data");

    data->length = json_array_size(jdata_data);
    data->ts_length = nav->observed_length;

    //n_obs
    if(opts->n_obs >= 0){
        data->n_obs = (opts->n_obs < data->length) ? opts->n_obs : data->length;
    } else {
        data->n_obs = data->length;
    }

    ssm_row_t **rows = malloc(data->length * sizeof (ssm_row_t *));
    if (rows==NULL) {
        ssm_print_err("Allocation impossible for ssm_data_row_t **");
        exit(EXIT_FAILURE);
    }

    data->length_nonan = 0;
    data->n_obs_nonan = 0;
    data->ind_nonan = ssm_u1_new(data->length);

    for (i=0; i< data->length; i++){
        rows[i] = malloc(sizeof (ssm_row_t));
        if (rows[i] == NULL) {
            ssm_print_err("Allocation impossible for ssm_data_row_t *");
            exit(EXIT_FAILURE);
        }

        json_t *jrow = json_array_get(jdata_data, i);

        json_t *jdate = json_object_get(jrow, "date");
        if(json_is_string(jdate)) {
            rows[i]->date = strdup(json_string_value(jdate));
        } else {
            snprintf(str, SSM_STR_BUFFSIZE, "error: data[%d].date is not a string\n", i);
            ssm_print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jtime = json_object_get(jrow, "time");
        if(json_is_number(jtime)) {
            rows[i]->time = (unsigned int) json_integer_value(jtime);
        } else {
            snprintf(str, SSM_STR_BUFFSIZE, "error: data[%d].time is not an integer\n", i);
            ssm_print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jobserved = json_object_get(jrow, "observed");
        rows[i]->ts_nonan_length = json_array_size(jobserved);

        if(rows[i]->ts_nonan_length){
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
                    snprintf(str, SSM_STR_BUFFSIZE, "error: data[%d].observed[%d] is not an integer\n", i, j);
                    ssm_print_err(str);
                    exit(EXIT_FAILURE);
                }
            }
        }

        rows[i]->values = ssm_load_jd1_new(jrow, "values");

        json_t *jreset = json_object_get(jrow, "reset");
        rows[i]->states_reset_length = json_array_size(jreset);

        if(rows[i]->states_reset_length){
            rows[i]->states_reset = malloc(rows[i]->states_reset_length * sizeof (ssm_state_t *));
            if (rows[i]->states_reset == NULL) {
                ssm_print_err("Allocation impossible for ssm_data_row_t.states_reset");
                exit(EXIT_FAILURE);
            }
        }

        for(j=0; j<rows[i]->states_reset_length; j++){
            json_t *jreset_j = json_array_get(jreset, j);
            if(json_is_number(jreset_j)) {
                int id = json_integer_value(jreset_j);
                rows[i]->states_reset[j] = nav->states[id];
            } else {
                snprintf(str, SSM_STR_BUFFSIZE, "error: data[%d].reset[%d] is not an integer\n", i, j);
                ssm_print_err(str);
                exit(EXIT_FAILURE);
            }
        }

        if(rows[i]->ts_nonan_length){
            data->ind_nonan[data->length_nonan] = i;
            data->length_nonan += 1;
            if(i< data->n_obs){
                data->n_obs_nonan += 1;
            }
        }
    }

    data->rows = rows;


    ssm_data_adapt_to_simul(data, jdata, nav, opts);

    return data;
}



void _ssm_row_free(ssm_row_t *row)
{
    free(row->date);
    if(row->ts_nonan_length){
        free(row->observed);
    }
    free(row->values);
    if(row->states_reset_length){
        free(row->states_reset);
    }
    free(row);
}


void ssm_data_free(ssm_data_t *data)
{
    int i;

    free(data->date_t0);
    free(data->ind_nonan);

    for(i=0; i< data->length; i++){
        _ssm_row_free(data->rows[i]);
    }
    free(data->rows);

    free(data);
}



/**
 * in case of simulation extend data or reset n_obs to take into
 * account opts->end ISO 8601 date
 */
void ssm_data_adapt_to_simul(ssm_data_t *data, json_t *jdata, ssm_nav_t *nav, ssm_options_t *opts)
{
    int i, n;

    if(strcmp("", opts->end)==0){
        return;
    }

    struct tm tm_start_date_t0;
    memset(&tm_start_date_t0, 0, sizeof(struct tm));

    strptime(data->date_t0, "%Y-%m-%d", &tm_start_date_t0);
    time_t t_start_date_t0 = mktime(&tm_start_date_t0);

    struct tm tm_end;
    memset(&tm_end, 0, sizeof(struct tm));
    strptime(opts->end, "%Y-%m-%d", &tm_end);
    time_t t_end = mktime(&tm_end);

    double delta_date_t0 = difftime(t_end, t_start_date_t0)/(24.0*60.0*60.0);
    if(delta_date_t0 < 0.0){
	ssm_print_err("end date is before t0");
	exit(EXIT_FAILURE);
    }

    unsigned int time_start;
    struct tm tm_start;
    memset(&tm_start, 0, sizeof(struct tm));
    time_t t_start;
    if(data->length){
        strptime(data->rows[data->length-1]->date, "%Y-%m-%d", &tm_start);
        time_start = data->rows[data->length-1]->time;
	t_start = mktime(&tm_start);
    } else {
        t_start = t_start_date_t0;
    }

    double delta = difftime(t_end, t_start)/(24.0*60.0*60.0);
    if(delta < 0.0){

	//there are data but t_end before the last data point. In this case we just adjuts n_obs and n_obs_nonan
	if(data->length){
	    data->n_obs = 0;
	    data->n_obs_nonan = 0;

	    while(data->rows[data->n_obs]->time <= delta_date_t0){
		if(data->rows[data->n_obs]->ts_nonan_length){
		    data->n_obs_nonan += 1;
		}
		data->n_obs += 1;
	    }
	    
	    //safety
	    data->n_obs = GSL_MIN(data->n_obs, data->length);
	    data->n_obs_nonan = GSL_MIN(data->n_obs, data->n_obs_nonan);

	    return;
	}

        ssm_print_err("end date is before t0");
        exit(EXIT_FAILURE);
    }

    int n_extra = (int) ceil(delta / (double) opts->freq);
    if(n_extra){
        int offset = data->length;
        data->length +=  n_extra;
        data->n_obs += n_extra;

        ssm_row_t **rows = realloc(data->rows, data->length * sizeof (ssm_row_t *));
        if (rows!=NULL) {
            data->rows = rows;
        } else {
            ssm_print_err("could not re-allocate memory for ssm_data_t rows");
            exit(EXIT_FAILURE);
        }

        time_t t = t_start;
        char iso_8601[] = "YYYY-MM-DD";
        double one_day_in_sec = 24.0*60.0*60.0;
        int inc = opts->freq * 24*60*60;
	json_t *jreset_all = json_object_get(jdata, "reset_all");
	int states_reset_length = json_array_size(jreset_all);

        for(n=0; n<n_extra; n++){
            ssm_row_t *row = malloc(sizeof (ssm_row_t));
            if (row == NULL) {
                ssm_print_err("Allocation impossible for ssm_data_row_t *");
                exit(EXIT_FAILURE);
            }

            t += inc;
            struct tm *tm;
            tm = gmtime(&t);

            strftime(iso_8601, sizeof(iso_8601), "%Y-%m-%d", tm);
            row->date = strdup(iso_8601);
            row->time = time_start + (unsigned int) difftime(t, t_start)/one_day_in_sec;
            row->ts_nonan_length = 0;


            row->states_reset_length = states_reset_length;
            row->states_reset = malloc(states_reset_length * sizeof (ssm_state_t *));
            if (row->states_reset == NULL) {
                ssm_print_err("Allocation impossible for ssm_data_row_t.states_reset");
                exit(EXIT_FAILURE);
            }
	    for(i=0; i<row->states_reset_length; i++){
		json_t *jreset_all_i = json_array_get(jreset_all, i);
		if(json_is_number(jreset_all_i)) {
		    int id = json_integer_value(jreset_all_i);
		    row->states_reset[i] = nav->states[id];
		} else {
		    char str[SSM_STR_BUFFSIZE];
		    snprintf(str, SSM_STR_BUFFSIZE, "error: reset_all[%d] is not an integer\n", i);
		    ssm_print_err(str);
		    exit(EXIT_FAILURE);
		}
	    }
	    
            row->values = NULL;
            data->rows[offset+n] = row;
        }
    }
}



ssm_calc_t *ssm_calc_new(json_t *jdata, ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, ssm_options_t *opts, int thread_id)
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

    unsigned long int seed;
    if(opts->flag_seed_time){
        seed = (unsigned) time(NULL);
    } else{
        seed=2;
    }
    calc->seed =  seed + opts->id; /*we ensure uniqueness of seed in case of parrallel runs*/

    calc->randgsl = gsl_rng_alloc(Type);
    gsl_rng_set(calc->randgsl, calc->seed + thread_id);

    /*******************/
    /* implementations */
    /*******************/
    int dim = _ssm_dim_X(nav);

    if (nav->implementation == SSM_ODE || nav->implementation == SSM_EKF){

        calc->T = gsl_odeiv2_step_rkf45;
        calc->control = gsl_odeiv2_control_y_new(opts->eps_abs, opts->eps_rel);
        calc->step = gsl_odeiv2_step_alloc(calc->T, dim);
        calc->evolve = gsl_odeiv2_evolve_alloc(dim);
        (calc->sys).function =  (nav->implementation == SSM_ODE) ? &ssm_step_ode: &ssm_step_ekf;
        (calc->sys).jacobian = NULL;
        (calc->sys).dimension= dim;
        (calc->sys).params= calc;

        if(nav->implementation == SSM_EKF){

            int can_run;

            if ( (nav->noises_off & (SSM_NO_DEM_STO)) && (nav->noises_off & (SSM_NO_WHITE_NOISE)) )  {
                calc->eval_Q = &ssm_eval_Q_no_dem_sto_no_env_sto;
                can_run = 0;
            } else if ((nav->noises_off & SSM_NO_DEM_STO) && !(nav->noises_off & SSM_NO_WHITE_NOISE)) {
                calc->eval_Q = &ssm_eval_Q_no_dem_sto;
                can_run = nav->par_noise->length;
            } else if (!(nav->noises_off & SSM_NO_DEM_STO) && (nav->noises_off & SSM_NO_WHITE_NOISE)) {
                calc->eval_Q = &ssm_eval_Q_no_env_sto;
                can_run = 1;
            } else {
                calc->eval_Q = &ssm_eval_Q_full;
                can_run = 1;
            }

            if(!(nav->noises_off & SSM_NO_DIFF)){
                can_run += nav->states_diff->length;
            }

            if(!can_run){
                ssm_print_err("Kalman methods must be used with at least one source of stochasticity in the process.");
                exit(EXIT_FAILURE);
            }

            int n_s = nav->states_sv_inc->length + nav->states_diff->length;
            int n_o = nav->observed_length;
            calc->_pred_error = gsl_vector_calloc(n_o);
            calc->_zero = gsl_vector_calloc(n_o);
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
            calc->_eval = gsl_vector_calloc(n_s);
            calc->_evec = gsl_matrix_calloc(n_s, n_s);
            calc->_w_eigen_vv = gsl_eigen_symmv_alloc(n_s);
        }

    } else if (nav->implementation == SSM_SDE){
        calc->y_pred = ssm_d1_new(dim);
    } else if (nav->implementation == SSM_PSR){
        ssm_psr_new(calc);
    }

    /**************************/
    /* multi-threaded sorting */
    /**************************/

    calc->J = fitness->J;
    calc->to_be_sorted = ssm_d1_new(fitness->J);
    calc->index_sorted = ssm_st1_new(fitness->J);

    /**************/
    /* references */
    /**************/
    calc->_par = NULL;
    calc->_nav = nav;

    /**************/
    /* covariates */
    /**************/

    json_t *jcovariates = json_object_get(jdata, "covariates");
    calc->covariates_length = json_array_size(jcovariates);

    if(calc->covariates_length){

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

        double freeze_forcing; // the time (in days) to freeze (i.e only take metadata from this time) (ignored if freeze_forcing < 0.0)
        double t_max; //t_max the highest possible time in days when interpolated metadata will be requested (negative values default to last point of metadata).

        //assess freeze_forcing and t_max...
        struct tm tm_start;
        memset(&tm_start, 0, sizeof(struct tm));
        strptime(data->date_t0, "%Y-%m-%d", &tm_start);
        time_t t_start = mktime(&tm_start);

        if(strcmp("", opts->end)!=0){
            struct tm tm_freeze;
            memset(&tm_freeze, 0, sizeof(struct tm));
            strptime(opts->freeze_forcing, "%Y-%m-%d", &tm_freeze);
            time_t t_freeze = mktime(&tm_freeze);
            freeze_forcing = difftime(t_freeze, t_start)/(24.0*60.0*60.0);
        } else {
            freeze_forcing = -1.0;
        }

        if(strcmp("", opts->end)!=0){
            struct tm tm_end;
            memset(&tm_end, 0, sizeof(struct tm));
            strptime(opts->end, "%Y-%m-%d", &tm_end);
            time_t t_end = mktime(&tm_end);
            t_max = difftime(t_end, t_start)/(24.0*60.0*60.0);
        } else {
            t_max = -1.0;
        }

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
                x_all[1] = GSL_MAX( GSL_MAX( t_max, ((data->n_obs>=1) ? (double) data->rows[data->n_obs-1]->time: 0.0) ),  x[size-1]);

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
    }
    return calc;
}


void ssm_calc_free(ssm_calc_t *calc, ssm_nav_t *nav)
{
    gsl_rng_free(calc->randgsl);

    if (nav->implementation == SSM_ODE  || nav->implementation == SSM_EKF){

        gsl_odeiv2_step_free(calc->step);
        gsl_odeiv2_evolve_free(calc->evolve);
        gsl_odeiv2_control_free(calc->control);

        if(nav->implementation == SSM_EKF){
            gsl_vector_free(calc->_pred_error);
            gsl_matrix_free(calc->_St);
            gsl_matrix_free(calc->_Stm1);
            gsl_matrix_free(calc->_Rt);
            gsl_matrix_free(calc->_Ht);
            gsl_matrix_free(calc->_Kt);
            gsl_matrix_free(calc->_Tmp_N_SV_N_TS);
            gsl_matrix_free(calc->_Tmp_N_TS_N_SV);
            gsl_matrix_free(calc->_Q);
            gsl_matrix_free(calc->_FtCt);
            gsl_matrix_free(calc->_Ft);
            gsl_vector_free(calc->_eval);
            gsl_matrix_free(calc->_evec);
            gsl_eigen_symmv_free(calc->_w_eigen_vv);
        }

    } else if (nav->implementation == SSM_SDE){
        free(calc->y_pred);
    } else if (nav->implementation == SSM_PSR){
        ssm_psr_free(calc);
    }

    free(calc->to_be_sorted);
    free(calc->index_sorted);

    if(calc->covariates_length){
        int k;
        for(k=0; k< calc->covariates_length; k++) {
            if(calc->spline[k]){
                gsl_spline_free(calc->spline[k]);
            }
            if(calc->acc[k]){
                gsl_interp_accel_free(calc->acc[k]);
            }
        }

        free(calc->spline);
        free(calc->acc);
    }

    free(calc);
}


ssm_calc_t **ssm_N_calc_new(json_t *jdata, ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, ssm_options_t *opts)
{
    int i;
    int n_threads = ssm_sanitize_n_threads(opts->n_thread, fitness);

    ssm_calc_t **calc = malloc(n_threads * sizeof (ssm_calc_t *));
    if (calc==NULL) {
        ssm_print_err("Allocation impossible for ssm_calc_t **");
        exit(EXIT_FAILURE);
    }

    for (i=0; i< n_threads; i++) {
        calc[i] = ssm_calc_new(jdata, nav, data, fitness, opts, i);
    }

    return calc;
}


void ssm_N_calc_free(ssm_calc_t **calc, ssm_nav_t *nav)
{
    int n;
    int threads_length = calc[0]->threads_length;

    for(n=0; n<threads_length; n++) {
        ssm_calc_free(calc[n], nav);
    }

    free(calc);
}



ssm_options_t *ssm_options_new(void)
{
    ssm_options_t *opts = malloc(sizeof(ssm_options_t));
    if (opts==NULL) {
        ssm_print_err("Allocation impossible for ssm_options_t *");
        exit(EXIT_FAILURE);
    }

    //alloc char *
    opts->freeze_forcing = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->root = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->next = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->interpolator = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->start = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->end = ssm_c1_new(SSM_STR_BUFFSIZE);
    opts->server = ssm_c1_new(SSM_STR_BUFFSIZE);

    //fill default
    opts->worker_algo = 0;

    opts->implementation = 0;
    opts->noises_off = 0;
    opts->print = 0;

    opts->id = 0;
    opts->flag_seed_time = 0;
    opts->flag_prior = 0;
    opts->dt = 0.25;
    opts->eps_abs = 1e-6;
    opts->eps_rel = 1e-3;
    strncpy(opts->freeze_forcing, "", SSM_STR_BUFFSIZE);
    strncpy(opts->root, ".", SSM_STR_BUFFSIZE);
    strncpy(opts->next, "", SSM_STR_BUFFSIZE);
    opts->n_thread = 1;
    opts->like_min = 1e-17;
    opts->J = 1;
    opts->n_obs = -1;
    strncpy(opts->interpolator, "linear", SSM_STR_BUFFSIZE);
    opts->n_obs = -1;
    opts->n_iter = 10;
    opts->a = 0.98;
    opts->b = 2;
    opts->L = 0.75;
    opts->m_switch = -1;
    opts->flag_ic_only = 0;
    opts->eps_switch = 50;
    opts->eps_max = 50.0;
    opts->flag_smooth = 0;
    opts->alpha = 0.02;
    opts->n_traj = 1000;
    opts->flag_tcp = 0;
    opts->flag_least_squares = 0;
    opts->size_stop = 1e-6;
    opts->freq = 1;
    strncpy(opts->start, "", SSM_STR_BUFFSIZE);
    strncpy(opts->end, "", SSM_STR_BUFFSIZE);
    strncpy(opts->server, "127.0.0.1", SSM_STR_BUFFSIZE);
    opts->flag_no_filter = 0;

    return opts;
}


void ssm_options_free(ssm_options_t *opts)
{
    free(opts->freeze_forcing);
    free(opts->root);
    free(opts->next);
    free(opts->interpolator);
    free(opts->start);
    free(opts->end);
    free(opts->server);

    free(opts);
}


ssm_fitness_t *ssm_fitness_new(ssm_data_t *data, ssm_options_t *opts)
{
    int i;

    ssm_fitness_t *fitness = malloc(sizeof(ssm_fitness_t));
    if (fitness==NULL) {
        ssm_print_err("Allocation impossible for ssm_fitness_t *");
        exit(EXIT_FAILURE);
    }

    fitness->J = opts->J;
    fitness->data_length = data->length;
    fitness->like_min = opts->like_min;
    fitness->log_like_min = log(fitness->like_min);

    fitness->ess_n = 0.0;
    fitness->log_like_n = 0.0;
    fitness->log_like = 0.0;

    fitness->weights = ssm_d1_new(fitness->J);
    fitness->select = ssm_u2_new(fitness->data_length, fitness->J);

    fitness->cum_status = malloc(fitness->J * sizeof (ssm_err_code_t));
    if(fitness->cum_status == NULL) {
        ssm_print_err("Allocation impossible for fitness->cum_status");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<fitness->J; i++){
        fitness->cum_status[i] = SSM_SUCCESS;
    }

    fitness->n_all_fail = 0;

    fitness->log_like_prev = 0.0;
    fitness->log_prior = 0.0;
    fitness->log_prior_prev = 0.0;

    return fitness;
}


void ssm_fitness_free(ssm_fitness_t *fitness)
{
    free(fitness->weights);
    ssm_u2_free(fitness->select, fitness->data_length);

    free(fitness->cum_status);

    free(fitness);
}


int _ssm_dim_X(ssm_nav_t *nav)
{
    int dim = nav->states_sv_inc->length + nav->states_diff->length;
    if(nav->implementation == SSM_EKF){
        dim += pow(dim, 2);
    }
    return dim;
}

ssm_X_t *ssm_X_new(ssm_nav_t *nav, ssm_options_t *opts)
{
    ssm_X_t *X = malloc(sizeof (ssm_X_t));
    if (X==NULL) {
        ssm_print_err("Allocation impossible for ssm_X_t");
        exit(EXIT_FAILURE);
    }

    X->length = _ssm_dim_X(nav);

    X->dt = 1.0/ ((double) round(1.0/opts->dt)); //IMPORTANT: for non adaptive time step methods, we ensure an integer multiple of dt in between 2 data points
    X->dt0 = X->dt;

    X->proj = ssm_d1_new(X->length);

    return X;
}

void ssm_X_free(ssm_X_t *X)
{
    free(X->proj);
    free(X);
}


ssm_X_t **ssm_J_X_new(ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_options_t *opts)
{
    int i;
    ssm_X_t **X = malloc(fitness->J * sizeof (ssm_X_t *));
    if (X==NULL) {
        ssm_print_err("Allocation impossible for ssm_X_t *");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<fitness->J; i++){
        X[i] = ssm_X_new(nav, opts);
    }

    return X;
}

void ssm_J_X_free(ssm_X_t **X, ssm_fitness_t *fitness)
{
    int i;

    for(i=0; i<fitness->J; i++){
        ssm_X_free(X[i]);
    }

    free(X);
}


ssm_X_t **ssm_D_X_new(ssm_data_t *data, ssm_nav_t *nav, ssm_options_t *opts)
{
    int i;
    ssm_X_t **X = malloc((data->length+1) * sizeof (ssm_X_t *));
    if (X==NULL) {
        ssm_print_err("Allocation impossible for ssm_X_t *");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<data->length+1; i++){
        X[i] = ssm_X_new(nav, opts);
    }

    return X;
}

void ssm_D_X_free(ssm_X_t **X, ssm_data_t *data)
{
    int i;

    for(i=0; i<data->length+1; i++){
        ssm_X_free(X[i]);
    }

    free(X);
}




ssm_X_t ***ssm_D_J_X_new(ssm_data_t *data, ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_options_t *opts)
{
    int i;
    ssm_X_t ***X = malloc((data->length+1) * sizeof (ssm_X_t **));
    if (X==NULL) {
        ssm_print_err("Allocation impossible for ssm_X_t **");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<data->length+1; i++){
        X[i] = ssm_J_X_new(fitness, nav, opts);
    }

    return X;
}

void ssm_D_J_X_free(ssm_X_t ***X, ssm_data_t *data, ssm_fitness_t *fitness)
{
    int i;

    for(i=0; i<data->length+1; i++){
        ssm_J_X_free(X[i], fitness);
    }

    free(X);
}


ssm_hat_t *ssm_hat_new(ssm_nav_t *nav)
{
    ssm_hat_t *hat = malloc(sizeof (ssm_hat_t));
    if (hat==NULL) {
        ssm_print_err("Allocation impossible for ssm_hat_t");
        exit(EXIT_FAILURE);
    }

    hat->states_length = nav->states_sv_inc->length + nav->states_diff->length;
    hat->states = ssm_d1_new(hat->states_length);
    hat->states_95 = ssm_d2_new(hat->states_length, 2);

    hat->remainders_length = nav->states_remainders->length;
    hat->remainders = ssm_d1_new(hat->remainders_length);
    hat->remainders_95 = ssm_d2_new(hat->remainders_length, 2);

    hat->observed_length = nav->observed_length;
    hat->observed = ssm_d1_new(hat->observed_length);
    hat->observed_95 = ssm_d2_new(hat->observed_length, 2);

    return hat;
}

void ssm_hat_free(ssm_hat_t *hat)
{
    free(hat->states);
    ssm_d2_free(hat->states_95, hat->states_length);

    free(hat->remainders);
    ssm_d2_free(hat->remainders_95, hat->remainders_length);

    free(hat->observed);
    ssm_d2_free(hat->observed_95, hat->observed_length);

    free(hat);
}


ssm_hat_t **ssm_D_hat_new(ssm_data_t *data, ssm_nav_t *nav)
{
    int i;

    ssm_hat_t **hat = malloc((data->length+1) * sizeof (ssm_hat_t *));
    if (hat==NULL) {
        ssm_print_err("Allocation impossible for ssm_hat_t *");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<data->length+1; i++){
        hat[i] = ssm_hat_new(nav);
    }

    return hat;
}


void ssm_D_hat_free(ssm_hat_t **hat, ssm_data_t *data)
{
    int i;

    for(i=0; i<data->length +1; i++){
        ssm_hat_free(hat[i]);
    }

    free(hat);
}


ssm_adapt_t *ssm_adapt_new(ssm_nav_t *nav, ssm_options_t * opts)
{
    ssm_adapt_t *a = malloc(sizeof (ssm_adapt_t));
    if(a == NULL) {
        ssm_print_err("allocation impossible for ssm_adapt_t");
        exit(EXIT_FAILURE);
    }

    a->ar = 1.0;
    a->ar_smoothed = 1.0;

    a->eps = 1.0;
    a->eps_max = opts->eps_max;
    a->eps_switch = opts->eps_switch;
    a->eps_a = opts->a;

    int min_switch = 5.0*pow(nav->theta_all->length, 2);
    if (opts->m_switch < 0) {
        a->m_switch = min_switch;
    } else {
        a->m_switch = opts->m_switch;
        if ( (a->m_switch < min_switch) && (nav->print & SSM_PRINT_WARNING)) {
            char str[SSM_STR_BUFFSIZE];
            snprintf(str, SSM_STR_BUFFSIZE, "warning: covariance switching iteration (%i) is smaller than proposed one (%i)", a->m_switch, min_switch);
            ssm_print_warning(str);
        }
    }

    a->flag_smooth = opts->flag_smooth;
    a->alpha = opts->alpha;

    a->mean_sampling = ssm_d1_new(nav->theta_all->length);
    a->var_sampling = gsl_matrix_calloc(nav->theta_all->length, nav->theta_all->length);

    return a;
}

void ssm_adapt_free(ssm_adapt_t *adapt)
{
    free(adapt->mean_sampling);
    gsl_matrix_free(adapt->var_sampling);

    free(adapt);
}
