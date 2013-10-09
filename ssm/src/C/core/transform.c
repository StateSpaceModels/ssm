/**************************************************************************
 *    This file is part of ssm.
 *
 *    ssm is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    ssm is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with ssm.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "ssm.h"

double ssm_f_id(double x)
{
    return x;
}


double ssm_f_log(double x)
{
    double safe = ( x > SSM_ZERO_LOG ) ? x : SSM_ZERO_LOG;
    return log(safe);
}

double ssm_f_inv_log(double x)
{
    return exp(x);
}

double ssm_f_logit(double x)
{
    //sanatize
    double safe = ( x > SSM_ZERO_LOG ) ? x : SSM_ZERO_LOG;
    safe = (safe < SSM_ONE_LOGIT ) ? safe : SSM_ONE_LOGIT;

    return log(safe/(1.0-safe));
}


double ssm_f_inv_logit(double x)
{
    if (x > 0) {
        return (1.0/(1.0+exp(-x)));
    } else {
        return (exp(x)/(1.0+exp(x)));
    }
}



double ssm_f_logit_ab(double x, double a, double b)
{
    if (a == b)
        return x; // nothing will happen in the transformed space for x, so no need to transform it
    else{
        double ratio = (x-a)/(b-x);
        if(ratio < SSM_ZERO_LOG){
            ratio = SSM_ZERO_LOG;
        } else if(ratio > (1.0/SSM_ZERO_LOG)) {
            ratio = 1.0/SSM_ZERO_LOG;
        }
        return log(ratio);
    }
}


double ssm_f_inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        if (x < 0) {
            return (b*exp(x)+a)/(1.0+exp(x));
        } else {
            return (b+a*exp(-x))/(1.0+exp(-x));
        };
    }
}

/**
 * derivative of ssm_f_id
 */
double ssm_f_der_id(double x)
{
    return 0;
}


/**
 * derivative of ssm_f_log
 */
double ssm_f_der_log(double x)
{
    return 1.0/x;
}



/**
 * derivative of ssm_f_inv_log
 */
double ssm_f_der_inv_log(double x)
{
    return exp(x);
}


/**
 * 2nd derivative of ssm_f_inv_log
 */
double ssm_f_der2_inv_log(double x)
{
    return exp(x);
}


/**
 * derivative of ssm_f_logit
 */
double ssm_f_der_logit(double x)
{
    return 1.0/(x-x*x);
}

/**
 * derivative of ssm_f_inv_logit
 */
double ssm_f_der_inv_logit(double x)
{
    if (x > 0) {
        return exp(-x)/pow(1.0 + exp(-x), 2.0);
    } else {
        return exp(x)/pow(1.0 + exp(x), 2.0);
    }
}


/**
 * 2nd derivative of ssm_f_inv_logit
 */
double ssm_f_der2_inv_logit(double x)
{
    if (x > 0) {
        return exp(-x)*(exp(-x)-1)/pow(1.0 + exp(-x), 3.0);
    } else {
        return exp(x)*(1-exp(x))/pow(1.0 + exp(x), 3.0);
    }
}


/**
 * derivative of ssm_f_logit_ab
 */
double ssm_f_der_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        return (b-a)/((x-a)*(b-x));
    }
}

/**
 * derivative of ssm_f_inv_logit_ab
 */
double ssm_f_der_inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        if (x > 0) {
            return (b-a)*exp(-x)/pow(exp(-x) + 1.0, 2.0);
        } else {
            return (b-a)*exp(x)/pow(exp(x) + 1.0, 2.0);
        }
    }
}


/**
 * 2nd derivative of ssm_f_inv_logit_ab
 */
double ssm_f_der2_inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        if (x > 0) {
            return (b-a)*exp(-x)*(exp(-x)-1)/pow(1.0 + exp(-x), 3.0);
        } else {
            return (b-a)*exp(x)*(1-exp(x))/pow(1.0 + exp(x), 3.0);
        }
    }
}


double ssm_f_user_par_id(double x, ssm_input_t *par, ssm_calc_t *calc)
{
    return x;
}

double ssm_f_2prior_id(double x, ssm_hat_t *hat, ssm_par_t *par, ssm_calc_t *calc, double t)
{
    return x;
}
