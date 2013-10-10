#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static ssm_nav_t *nav;
static ssm_options_t *opts;
static ssm_input_t *input;
static ssm_var_t *var;
static ssm_X_t *X;

void test_inputs__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("datapackage.json"));
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);

    input = ssm_input_new(jparameters, nav);
    var = ssm_var_new(jparameters, nav);
    X = ssm_X_new(nav, opts);

}

void test_inputs__cleanup(void)
{
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);

    ssm_input_free(input);
    ssm_var_free(var);
    ssm_X_free(X);
}

void test_inputs__input_new(void)
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


void test_inputs__var_new(void)
{
    int i, j;
    double expected[10][10] = {
        {0.01, 0,    0,    0,    0,    0,    0,    0,    0,    0},
        {0,    0.03, 0,    0,    0,    0,    0,    0,    0,    0},
        {0,    0,    0.02, 0,    0,    0,    0,    0,    0,    0},
        {0,    0,    0,    0.04, 0,    0.01, 0,    0,    0,    0},
        {0,    0,    0,    0,    0.02, 0,    0,    0,    0,    0},
        {0,    0,    0,    0.01, 0,    0.02, 0,    0,    0,    0},
        {0,    0,    0,    0,    0,    0,    0.02, 0,    0,    0},
        {0,    0,    0,    0,    0,    0,    0,    0.02, 0,    0},
        {0,    0,    0,    0,    0,    0,    0,    0,    0.02, 0},
        {0,    0,    0,    0,    0,    0,    0,    0,    0,    0.03}
    };

    for(i=0; i<nav->theta_all->length; i++){
        for(j=0; j<nav->theta_all->length; j++){
            cl_check(gsl_matrix_get(var, nav->theta_all->p[i]->offset_theta, nav->theta_all->p[j]->offset_theta) == expected[i][j]);
        }
    }

}

void test_inputs__X_new(void)
{
    int i;
    cl_check(X->length == 4+2+2);
    cl_check(X->dt0 == 0.25);
    cl_check(X->dt == X->dt0);

    for(i=0; i<X->length; i++){
        cl_check(X->proj[i] == 0.0);
    }
}
