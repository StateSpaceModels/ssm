#include "clar.h"
#include <ssm.h>

static json_t *jparameters;
static ssm_nav_t *nav;
static ssm_options_t *opts;

void test_nav__initialize(void)
{
    jparameters = ssm_load_json_file(cl_fixture("theta.json"));    
    opts = ssm_options_new();
    nav = ssm_nav_new(jparameters, opts);
}

void test_nav__cleanup(void)
{
    json_decref(jparameters);
    ssm_options_free(opts);
    ssm_nav_free(nav);
}


void test_nav__nav_it_theta(void)
{        
    int i;

    cl_check(nav->theta_all->length == 10);
    cl_check(nav->theta_no_icsv_no_icdiff->length == 5);
    cl_check(nav->theta_icsv_icdiff->length == 5);

    char *expected_names[] = {"pr_I_nyc", "pr_I_paris", "pr_S_paris", "r0_nyc", "r0_paris", "pr_v", "rep_all_CDC_inc", "rep_all_google_inc", "rep_nyc_CDC_inc", "rep_paris_CDC_prev"};
    for(i=0; i<nav->theta_all->length; i++){
	cl_assert_equal_s(nav->theta_all->p[i]->name, expected_names[i]);
    }

    char *expected_names_ic[] = {"pr_I_nyc", "pr_I_paris", "pr_S_paris", "r0_nyc", "r0_paris"};
    for(i=0; i<nav->theta_icsv_icdiff->length; i++){
	cl_assert_equal_s(nav->theta_icsv_icdiff->p[i]->name, expected_names_ic[i]);
    }

    char *expected_names_no[] = {"pr_v", "rep_all_CDC_inc", "rep_all_google_inc", "rep_nyc_CDC_inc", "rep_paris_CDC_prev"};
    for(i=0; i<nav->theta_no_icsv_no_icdiff->length; i++){
	cl_assert_equal_s(nav->theta_no_icsv_no_icdiff->p[i]->name, expected_names_no[i]);
    }

}


void test_nav__nav_it_theta_opts_no_diff(void)
{        
    int i;

    ssm_options_t *opts_nd = ssm_options_new();
    opts_nd->noises_off = SSM_NO_DIFF;
    ssm_nav_t *nav_nd = ssm_nav_new(jparameters, opts_nd);

    cl_check(nav_nd->theta_no_icsv_no_icdiff->length == 7);
    cl_check(nav_nd->theta_icsv_icdiff->length == 3);

    char *expected_names[] = {"pr_I_nyc", "pr_I_paris", "pr_S_paris"};
    for(i=0; i<nav_nd->theta_icsv_icdiff->length; i++){
	cl_assert_equal_s(nav_nd->theta_icsv_icdiff->p[i]->name, expected_names[i]);
    }

    char *expected_names_no[] = {"r0_nyc", "r0_paris", "pr_v", "rep_all_CDC_inc", "rep_all_google_inc", "rep_nyc_CDC_inc", "rep_paris_CDC_prev"};
    for(i=0; i<nav_nd->theta_no_icsv_no_icdiff->length; i++){       
	cl_assert_equal_s(nav_nd->theta_no_icsv_no_icdiff->p[i]->name, expected_names_no[i]);
    }

    ssm_options_free(opts_nd);
    ssm_nav_free(nav_nd);
}
