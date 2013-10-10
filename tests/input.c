#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static ssm_nav_t *nav;
static ssm_options_t *opts;
static ssm_input_t *input;

void test_input__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("datapackage.json"));    
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);

    input = ssm_input_new(jparameters, nav);
}

void test_input__cleanup(void)
{
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);

    ssm_input_free(input);
}


void test_input__input_new(void)
{        
    int i;

    double expected[] = {
	1e-05, //I_nyc
	1e-05, //I_paris
	0.07,  //S_nyc
	0,     //S_paris
	0.1,   //sto
	20,    //r0_nyc
	20,    //r0_paris
	11,    //v
	0,     //vol
	0.1,   //phi
	0.6,   //rep_all_CDC_inc
	0.6,   //rep_all_google_inc
	0.6,   //rep_nyc_CDC_inc
	0.6    //rep_paris_CDC_prev
    };

    for(i=0; i<nav->par_all->length; i++){
	cl_check(gsl_vector_get(input, nav->par_all->p[i]->offset) == expected[i]);
    }
}
