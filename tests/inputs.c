#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static json_t *jdata;
static ssm_options_t *opts;
static ssm_nav_t *nav;
static ssm_data_t *data;
static ssm_calc_t *calc;
static ssm_fitness_t *fitness;
static ssm_input_t *input;
static ssm_var_t *var;
static ssm_X_t *X;
static ssm_par_t *par;
static ssm_theta_t *theta;

void test_inputs__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("datapackage.json"));
    jdata = ssm_load_json_file(cl_fixture("data.json"));
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);
    data = ssm_data_new(jdata, nav, opts);
    fitness = ssm_fitness_new(data, opts);
    calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);

    input = ssm_input_new(jparameters, nav);
    var = ssm_var_new(jparameters, nav);
    X = ssm_X_new(nav, opts);
    par = ssm_par_new(input, calc, nav);
    theta = ssm_theta_new(input, nav);
}

void test_inputs__cleanup(void)
{
    ssm_calc_free(calc, nav);
    json_decref(jdata);
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);
    ssm_data_free(data);
    ssm_fitness_free(fitness);

    ssm_input_free(input);
    ssm_var_free(var);
    ssm_X_free(X);
    ssm_par_free(par);
    ssm_theta_free(theta);
}

void test_inputs__input_new(void)
{
    int i;
    double expected[] = {
        5,     //I_nyc
        1e-05, //I_paris
        0,     //S_nyc (follow are initialized at 0)
        0.07,  //S_paris
        0.1,   //sto
        20,    //r0_nyc
        20,    //r0_paris
        11,    //v
        0.1,   //vol
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

void test_inputs__par_new(void)
{
    int i;
    double expected[] = {
        pow(10,-5)*1000000, //I_nyc
        1e-05*1000000,      //I_paris
        0.07*1000000,       //S_nyc
        0.07*1000000,       //S_paris
        0.1,                //sto
        20,                 //r0_nyc
        20,                 //r0_paris
        1.0/11.0,           //v
        0.1,                //vol
        0.1,                //phi
        0.6,                //rep_all_CDC_inc
        0.6,                //rep_all_google_inc
        0.6,                //rep_nyc_CDC_inc
        0.6                 //rep_paris_CDC_prev
    };

    for(i=0; i<nav->par_all->length; i++){
        cl_check(gsl_vector_get(par, nav->par_all->p[i]->offset) == expected[i]);
    }
}


void test_inputs__theta_new(void)
{
    int i;
    double expected[] = {
        ssm_f_logit_ab(5, 4, 6),               //I_nyc
        ssm_f_logit_ab(1e-05, 1e-6, 1e-4),     //I_paris
        ssm_f_logit_ab(0.07, 0.04, 0.09),      //S_paris
        ssm_f_logit_ab(20.0, 15.0, 35.0),      //r0_nyc
        ssm_f_logit_ab(20.0, 15.0, 35.0),      //r0_paris
        ssm_f_log(11.0),                       //v
        ssm_f_logit_ab(0.6, 0.5, 0.8),         //rep_all_CDC_inc
        ssm_f_logit_ab(0.6, 0.5, 0.8),         //rep_all_google_inc
        ssm_f_logit_ab(0.6, 0.5, 0.8),         //rep_nyc_CDC_inc
        ssm_f_logit_ab(0.6, 0.5, 0.8)          //rep_paris_CDC_prev
    };

    for(i=0; i<nav->theta_all->length; i++){
	cl_check(gsl_vector_get(theta, nav->theta_all->p[i]->offset_theta) == expected[i]);
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

