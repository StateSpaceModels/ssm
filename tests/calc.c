#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static json_t *jdata;
static ssm_nav_t *nav;
static ssm_options_t *opts;
static ssm_data_t *data;
static ssm_fitness_t *fitness;
static ssm_calc_t *calc;

void test_calc__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("package.json"));
    jdata = ssm_load_json_file(cl_fixture(".data.json"));
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);
    data = ssm_data_new(jdata, nav, opts);
    fitness = ssm_fitness_new(data, opts);
    calc = ssm_calc_new(jdata, nav, data, fitness, opts, 0);
}

void test_calc__cleanup(void)
{
    ssm_calc_free(calc, nav);
    json_decref(jdata);
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);
    ssm_data_free(data);
    ssm_fitness_free(fitness);
}

void test_calc__calc_new(void)
{
    int j;

    cl_check(calc->seed == 2);
    cl_check(calc->covariates_length == 10);
    cl_check(calc->J == 1);

    for(j=0; j<calc->J; j++){
        cl_check(calc->to_be_sorted[j] == 0.0);
        cl_check(calc->index_sorted[j] == 0);
    }
}

void test_calc__covariates(void)
{
    int i;
    for(i=0; i<calc->covariates_length; i++){
        cl_assert_equal_s(gsl_spline_name(calc->spline[i]), "linear");
    }

    //par_forced: ['N_nyc', 'N_paris', 'mu_b_nyc', 'mu_b_paris', 'mu_d_nyc', 'mu_d_paris', 'prop_all_CDC_inc', 'prop_all_google_inc', 'prop_nyc_CDC_inc', 'prop_paris_CDC_prev']
    for(i=0; i<data->length; i++){
        cl_check(gsl_spline_eval(calc->spline[0], data->rows[i]->time, calc->acc[0]) == 1000000);
        cl_check(gsl_spline_eval(calc->spline[1], data->rows[i]->time, calc->acc[1]) == 1000000);

        cl_check(gsl_spline_eval(calc->spline[2], data->rows[i]->time, calc->acc[2]) == 0.00027);
        cl_check(gsl_spline_eval(calc->spline[3], data->rows[i]->time, calc->acc[3]) == 0.00027);
        cl_check(gsl_spline_eval(calc->spline[4], data->rows[i]->time, calc->acc[4]) == 0.00027);
        cl_check(gsl_spline_eval(calc->spline[5], data->rows[i]->time, calc->acc[5]) == 0.00027);

        cl_check(gsl_spline_eval(calc->spline[6], data->rows[2]->time, calc->acc[6]) == 1.0);
        cl_check(gsl_spline_eval(calc->spline[7], data->rows[3]->time, calc->acc[7]) == 1.0);
        cl_check(gsl_spline_eval(calc->spline[8], data->rows[2]->time, calc->acc[8]) == 1.0);
        cl_check(gsl_spline_eval(calc->spline[9], data->rows[3]->time, calc->acc[9]) == 1.0);
    }
}
