#include "clar.h"
#include <ssm.h>


static int parameters_length; 
static int states_length; 
static ssm_parameter_t **parameters;
static ssm_state_t **states;

void test_iterators__initialize(void)
{
    parameters = _ssm_parameters_new(&parameters_length);
    states = _ssm_states_new(&states_length, parameters);
}

void test_iterators__cleanup(void)
{
    int i;
    for(i=0; i<parameters_length; i++){
        _ssm_parameter_free(parameters[i]);
    }
    free(parameters);

    for(i=0; i<states_length; i++){
        _ssm_state_free(states[i]);
    }
    free(states);
}

void test_iterators__states_sv(void)
{    
    int i;
    ssm_it_states_t *states_sv = ssm_it_states_sv_new(states);
    char *expected_names[] = {"I_nyc", "I_paris", "S_nyc", "S_paris"};
    cl_check(states_sv->length == sizeof(expected_names)/sizeof(*expected_names));
    for(i=0; i<states_sv->length; i++){
	cl_assert_equal_s(states_sv->p[i]->name, expected_names[i]);
    }
    _ssm_it_states_free(states_sv);
}

void test_iterators__par_all_and_ssm_in_par(void)
{    
    int i;
    ssm_it_parameters_t *par_all = ssm_it_parameters_all_new(parameters);
    char *expected_names[] = {"pr_I_nyc", "pr_I_paris", "S_nyc", "pr_S_paris", "sto", "r0_nyc", "r0_paris", "pr_v", "vol", "phi", "rep_all_CDC_inc", "rep_all_google_inc", "rep_nyc_CDC_inc", "rep_paris_CDC_prev"};

    cl_check(par_all->length == sizeof(expected_names)/sizeof(*expected_names));
    for(i=0; i<par_all->length; i++){
	cl_assert_equal_s(par_all->p[i]->name, expected_names[i]);
    }

    cl_assert(ssm_in_par(par_all, "r0_paris") == 1);
    cl_assert(ssm_in_par(par_all, "x") == 0);

    _ssm_it_parameters_free(par_all);
}
