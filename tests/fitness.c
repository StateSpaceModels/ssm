#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static json_t *jdata;
static ssm_nav_t *nav;
static ssm_options_t *opts;
static ssm_data_t *data;
static ssm_fitness_t *fitness;

void test_fitness__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("datapackage.json"));
    jdata = ssm_load_json_file(cl_fixture("data.json"));
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);
    data = ssm_data_new(jdata, nav, opts);
    fitness = ssm_fitness_new(data, opts);
}

void test_fitness__cleanup(void)
{
    json_decref(jdata);
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);
    ssm_data_free(data);
    ssm_fitness_free(fitness);
}

void test_fitness__fitness_new(void)
{
    int i, j;

    cl_check(fitness->J == 1);
    cl_check(fitness->data_length == data->length);

    cl_check(fitness->like_min == 1e-17);
    cl_check(fitness->log_like_min == log(1e-17));

    cl_check(fitness->ess_n == 0.0);
    cl_check(fitness->log_like_n == 0.0);
    cl_check(fitness->log_like == 0.0);

    cl_check(fitness->n_all_fail == 0);
    cl_check(fitness->log_like_prev == 0.0);
    cl_check(fitness->log_prior == 0.0);
    cl_check(fitness->log_prior_prev == 0.0);

    for(j=0; j<fitness->J; j++){
        cl_check(fitness->weights[j] == 0.0);
        cl_check(fitness->cum_status[j] == SSM_SUCCESS);
    }

    for(i=0; i<fitness->data_length; i++){
        for(j=0; j<fitness->J; j++){
            cl_check(fitness->select[i][j] == 0);
        }
    }

}
