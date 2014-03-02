#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static json_t *jdata;
static ssm_nav_t *nav;
static ssm_options_t *opts;
static ssm_data_t *data;

void test_data__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("theta.json"));
    jdata = ssm_load_json_file(cl_fixture(".data.json"));
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);
    data = ssm_data_new(jdata, nav, opts);
}

void test_data__cleanup(void)
{
    json_decref(jdata);
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);
    ssm_data_free(data);
}

void test_data__data_new(void)
{
    cl_check(data->length == 52);
    cl_check(data->ts_length == 4);
    cl_check(data->n_obs == 52);
    cl_check(data->n_obs_nonan == 47);
    cl_check(data->length_nonan == 47);
    cl_assert_equal_s(data->date_t0, "2012-07-26");
}

void test_data__ind_nonan(void)
{
    int i, j;
    j = 0;
    cl_check(data->ind_nonan[0] == 0);
    cl_check(data->ind_nonan[1] == 1);
    for(i=2; i<47; i++){
        cl_check(data->ind_nonan[i] == i+5);
    }
}

void test_data__rows(void)
{
    ssm_row_t *row = data->rows[0];
    cl_assert_equal_s(row->date, "2012-08-02");
    cl_check(row->time == 7);
    cl_check(row->ts_nonan_length == 4);
    cl_check(row->states_reset_length == 2);

    cl_assert_equal_s(row->observed[0]->name, "all_CDC_inc");
    cl_assert_equal_s(row->observed[1]->name, "paris_CDC_prev");
    cl_assert_equal_s(row->observed[2]->name, "all_google_inc");
    cl_assert_equal_s(row->observed[3]->name, "nyc_CDC_inc");

    cl_check(row->values[0] == 6);
    cl_check(row->values[1] == 5);
    cl_check(row->values[2] == 10);
    cl_check(row->values[3] == 9);

    cl_assert_equal_s(row->states_reset[0]->name, "Inc_in_nyc");
    cl_assert_equal_s(row->states_reset[1]->name, "Inc_out");

    row = data->rows[29];
    cl_assert_equal_s(row->date, "2013-02-21");
    cl_check(row->time == 210);
    cl_check(row->ts_nonan_length == 3);
    cl_check(row->states_reset_length == 2);

    cl_assert_equal_s(row->observed[0]->name, "paris_CDC_prev");
    cl_assert_equal_s(row->observed[1]->name, "all_google_inc");
    cl_assert_equal_s(row->observed[2]->name, "nyc_CDC_inc");

    cl_check(row->values[0] == 854);
    cl_check(row->values[1] == 797);
    cl_check(row->values[2] == 2);

    cl_assert_equal_s(row->states_reset[0]->name, "Inc_in_nyc");
    cl_assert_equal_s(row->states_reset[1]->name, "Inc_out");

}


void test_data__ssm_data_adapt_to_simul(void)
{
    strncpy(opts->end, "2014-07-25", SSM_STR_BUFFSIZE);
    int prev_length = data->length;
    ssm_data_adapt_to_simul(data, jdata, nav, opts);

    cl_assert_equal_s(data->rows[prev_length]->date, "2013-07-26");
    cl_check(data->rows[prev_length]->time == 365);

    cl_assert_equal_s(data->rows[data->length-1]->date, "2014-07-25");
    cl_check(data->rows[data->length-1]->time == 729);
}


void test_data__ssm_data_adapt_to_simul_non_extra_domain(void)
{
    strncpy(opts->end, "2012-12-06", SSM_STR_BUFFSIZE);
    ssm_data_adapt_to_simul(data, jdata, nav, opts);

    cl_check(data->n_obs == 19);
    cl_check(data->n_obs_nonan == 14);
}
