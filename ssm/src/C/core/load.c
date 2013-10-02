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

json_t *load_json(void)
{
    char *buffer;
    buffer = malloc(SSM_BUFFER_SIZE* sizeof(char));

    fgets(buffer, SSM_BUFFER_SIZE, stdin);
    //    printf("%s\n", buffer);

    json_t *root;
    json_error_t error;

    root = json_loads(buffer, 0, &error);
    if(!root) {
        char str[SSM_STR_BUFFSIZE];
        sprintf(str, "could not parse parameters datapackage\nerror: on line %d: %s\n", error.line, error.text);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //  buffer = json_dumps(root, 0);

    free(buffer);

    return root;
}


/**
 *load constants defined as global variable
 */
json_t *load_settings(const char *path)
{
    json_error_t settings_error;
    json_t *settings = json_load_file(path, 0, &settings_error);
    if(!settings) {
        print_err(settings_error.text);
        exit(EXIT_FAILURE);
    }
    
    return settings;
}



void ssm_input2par(ssm_par_t *par, ssm_input_t *input, ssm_calc_t *calc, ssm_nav_t *nav)
{   
    ssm_it_parameters_t *it = nav->theta_all;
    
    int i;
    
    for(i=0; i< it->length; i++){
	gsl_vector_set(par, i, it->p[i]->f_user2par(gsl_vector_get(input, i), input, calc));
    }
}


void ssm_par2X(ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav)
{
    int i;

    ssm_it_states_t *sv = nav->states_sv;
    ssm_it_states_t *inc = nav->states_inc;
    ssm_it_states_t *diff = nav->states_diff;
    
    for(i=0; i<sv->length; i++){
	gsl_vector_set(X, sv->p[i]->offset, gsl_vector_get(par, sv->p[i]->ic->offset));
    }

    for(i=0; i<inc->length; i++){
	gsl_vector_set(X, inc->p[i]->offset, 0.0);
    }

    for(i=0; i<diff->length; i++){
	gsl_vector_set(X, diff->p[i]->offset, gsl_vector_get(par, diff->p[i]->ic->offset));
    }

}


unsigned int *ssm_load_ju1_new(json_t *container, char *name)
{
    char str[SSM_STR_BUFFSIZE];
    int i;
    json_t *array = json_object_get(container, name);
    unsigned int *tab = malloc(json_array_size(array) * sizeof (unsigned int));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);

        if(json_is_number(array_i)){
	    tab[i] = (unsigned int) json_integer_value(array_i);
        } else if(json_is_null(array_i)) {
	    tab[i] = NAN;
	} else {
            sprintf(str, "error: %s[%d] is not a number\n", name, i);
            print_err(str);
            exit(EXIT_FAILURE);
	}
    }

    return tab;
}


double *ssm_load_jd1_new(json_t *container, char *name)
{
    char str[SSM_STR_BUFFSIZE];
    int i;
    json_t *array = json_object_get(container, name);
    double *tab = malloc(json_array_size(array) * sizeof (double));
    if(tab==NULL){
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);

        if(json_is_number(array_i)) {
            tab[i] = json_number_value(array_i);
        } else if(json_is_null(array_i)) {
            tab[i] = NAN;
        } else {
            sprintf(str, "error: %s[%d] is not a number nor null\n", name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
    }

    return tab;
}


char **ssm_load_jc1_new(json_t *container, const char *name)
{

    char str[SSM_STR_BUFFSIZE];
    int i;
    json_t *array = json_object_get(container, name);
    char **tab = malloc(json_array_size(array) * sizeof (char *));
    if(tab==NULL){
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);

        if(json_is_string(array_i)) {
            tab[i] = strdup(json_string_value(array_i));
        } else {
            sprintf(str, "error: %s[%d] is not a string\n", name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }	
    }

    return tab;
}
