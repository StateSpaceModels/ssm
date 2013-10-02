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

ssm_var_t *ssm_var_new(ssm_nav_t *nav, json_t *parameters)
{
    gsl_matrix *m = gsl_matrix_calloc(nav->theta_all->length, nav->theta_all->length);


    int i,j, index;
    ssm_it_parameters_t *it = nav->theta_all;
    json_t *resource = json_object_get(parameters, "resource");

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
                            if(json_is_number(jcov_ijjval)) {
                                gsl_matrix_set(m, i, j, json_number_value(jcov_ijjval));
                            } else {
                                char str[SSM_STR_BUFFSIZE];
                                sprintf(str, "error: parameters.covariance.%s.%s is not a number\n", it->p[i]->name, it->p[j]->name);
                                print_err(str);
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
    int i;

    ssm_it_states_t *it = malloc(sizeof (ssm_it_states_t));
    if (it == NULL) {
        print_err("Allocation impossible for ssm_it_states_t *");
        exit(EXIT_FAILURE);
    }


    it->length = length;
    if(length){
        it->p = malloc(length * sizeof (ssm_state_t *));
        if (it->p == NULL) {
            print_err("Allocation impossible for ssm_it_states_t *");
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
    int i;

    ssm_it_parameters_t *it = malloc(sizeof (ssm_it_parameters_t));
    if (it == NULL) {
        print_err("Allocation impossible for ssm_it_parameters_t *");
        exit(EXIT_FAILURE);
    }

    it->length = length;
    if(length){
        it->p = malloc(length * sizeof (ssm_parameter_t *));
        if (it->p == NULL) {
            print_err("Allocation impossible for ssm_it_parameters_t *");
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


ssm_nav_t *ssm_nav_new(json_t jparameters, ...)
{
    ssm_nav_t *nav = malloc(sizeof (ssm_nav_t));
    if (nav == NULL) {
        print_err("Allocation impossible for ssm_nav_t *");
        exit(EXIT_FAILURE);
    }

    nav->parameters = ssm_parameters_new(&nav->parameters_length);
    nav->states = ssm_states_new(&nav->states_length, nav->parameters);
    nav->observed = ssm_observed_new(&nav->observed_length);

    nav->states_sv = ssm_it_states_sv_new(nav->states);
    nav->states_remainders = ssm_it_states_remainders_new(nav->states);
    nav->states_inc = ssm_it_incs_sv_new(nav->states);
    nav->states_diff = ssm_it_diff_sv_new(nav->states);

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

    int index;
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
                            print_err(str);
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
                                nav->theta_no_icsv_no_icdiff->p[nav->no_icsv_no_icdiff->length] = nav->par_all->p[i];
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


ssm_data_t *ssm_data_new(json_t *jdata, ssm_nav_t *nav, ssm_options *opts)
{
    char str[SSM_STR_BUFFSIZE];
    int i, j;

    ssm_data_t *p_data = malloc(sizeof (ssm_data_t));
    if (p_data==NULL) {
        print_err("Allocation impossible for ssm_data_t");
        exit(EXIT_FAILURE);
    }

    p_data->dates_t0 = ssm_load_jc1_new(jdata, "starts");

    json_t *jdata_data = json_object_get(jdata, "data");

    p_data->length = json_array_size(jdata_data);
    p_data->ts_length = nav->observed_length;

    ssm_data_row_t **rows = malloc(p_data->length * sizeof (ssm_data_row_t *));
    if (rows==NULL) {
        print_err("Allocation impossible for ssm_data_row_t **");
        exit(EXIT_FAILURE);
    }

    p_data->length_nonan = 0;
    p_data->ind_nonan = ssm_d1_new(p_data->length);

    for (i=0; i< p_data->length; i++){
        rows[i] = malloc(sizeof (p_data->length * (ssm_data_row_t)));
        if (rows[i] == NULL) {
            print_err("Allocation impossible for ssm_data_row_t *");
            exit(EXIT_FAILURE);
        }

        json_t *jrow = json_array_get(jdata_data, i);

        json_t *jdate = json_object_get(jrow, "date");
        if(json_is_string(jdate)) {
            rows[i]->date = strdup(json_string_value(jdate));
        } else {
            sprintf(str, "error: data[%d].date is not a string\n", i);
            print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jtime = json_object_get(jrow, "time");
        if(json_is_number(jtime)) {
            rows[i]->time = (unsigned int) json_integer_value(jtime);
        } else {
            sprintf(str, "error: data[%d].time is not an integer\n", i);
            print_err(str);
            exit(EXIT_FAILURE);
        }

        json_t *jobserved = json_object_get(jrow, "observed");
        rows[i]->ts_nonan_length = json_array_size(jobserved);
        rows[i]->observed = malloc(rows[i]->ts_nonan_length * sizeof (ssm_observed_t *));
        if (rows[i]->observed == NULL) {
            print_err("Allocation impossible for ssm_data_row_t.observed");
            exit(EXIT_FAILURE);
        }

        for(j=0; j<rows[i]->ts_nonan_length; j++){

            json_t *jobserved_j = json_array_get(jobserved, j);
            if(json_is_number(jobserved_j)) {
                int id = json_integer_value(jobserved_j);
                rows[i]->observed[j] = nav->observed[id];
            } else {
                sprintf(str, "error: data[%d].observed[%d] is not an integer\n", i, j);
                print_err(str);
                exit(EXIT_FAILURE);
            }
        }

        rows[i]->values = ssm_load_jd1_new(jrow, "values");

        json_t *jreset = json_object_get(jrow, "reset");
        rows[i]->states_reset_length = json_array_size(jreset);
        rows[i]->states_reset = malloc(rows[i]->states_reset_length * sizeof (ssm_state_t *));
        if (rows[i]->states_reset == NULL) {
            print_err("Allocation impossible for ssm_data_row_t.states_reset");
            exit(EXIT_FAILURE);
        }

        for(j=0; j<rows[i]->states_reset_length; j++){

            json_t *jreset_j = json_array_get(jreset, j);
            if(json_is_number(jreset_j)) {
                int id = json_integer_value(jreset_j);
                rows[i]->reset[j] = nav->states[id];
            } else {
                sprintf(str, "error: data[%d].reset[%d] is not an integer\n", i, j);
                print_err(str);
                exit(EXIT_FAILURE);
            }
        }

        if(rows[i]->ts_nonan_length){
            p_data->ind_nonan[p_data->length_nonan] = i;
            p_data->length_nonan += 1;
        }
    }

    p_data->rows = rows;

    //n_obs
    if(opts->n_obs >= 0){
        p_data->n_obs = (opts->n_obs < p_data->length) ? opts->n_obs : p_data->length;
    } else {
        p_data->n_obs = p_data->length;
    }

    return p_data;
}
