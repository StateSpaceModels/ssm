#include "clar.h"
#include <ssm.h>

static int observed_length; 
static ssm_observed_t **observed;

void test_observed__initialize(void)
{
    observed = _ssm_observed_new(&observed_length);
}

void test_observed__cleanup(void)
{
    int i;
    for(i=0; i<observed_length; i++){
        _ssm_observed_free(observed[i]);
    }
    free(observed);
}

void test_observed__observed_new(void)
{    
    int i;
    char *expected_names[] = {"all_CDC_inc", "paris_CDC_prev", "all_google_inc", "nyc_CDC_inc"};

    cl_check(observed_length == sizeof(expected_names)/sizeof(*expected_names));
    for(i=0; i<observed_length; i++){
	cl_assert(observed[i]->offset == i);
	cl_assert_equal_s(observed[i]->name, expected_names[i]);
    }
}
