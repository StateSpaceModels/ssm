#include "clar.h"
#include <ssm.h>

static int parameters_length; 
static ssm_parameter_t **parameters;

void test_parameters__initialize(void)
{
    parameters = _ssm_parameters_new(&parameters_length);
}

void test_parameters__cleanup(void)
{
    int i;
    for(i=0; i<parameters_length; i++){
        _ssm_parameter_free(parameters[i]);
    }
}

void test_parameters__parameters_new(void)
{    
    int i;
    char *expected_names[] = {"I_nyc", "I_paris", "S_nyc", "S_paris", "sto", "r0_nyc", "r0_paris", "v", "vol", "phi", "rep_all_CDC_inc", "rep_all_google_inc", "rep_nyc_CDC_inc", "rep_paris_CDC_prev"};

    cl_check(parameters_length == sizeof(expected_names)/sizeof(*expected_names));
    for(i=0; i<parameters_length; i++){
	cl_assert(parameters[i]->offset == i);
	cl_assert(parameters[i]->offset_theta == -1);
	cl_assert_equal_s(parameters[i]->name, expected_names[i]);
    }

}
