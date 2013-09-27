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

#ifndef SSM_H
#define SSM_H

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>

#include <getopt.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#include <gsl/gsl_spline.h>

#include <gsl/gsl_sort.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <jansson.h> //json

//parallel computing ability
#include <zmq.h>
#include <pthread.h>


enum ssm_implementations {SSM_ODE, SSM_SDE, SSM_PSR};
enum ssm_noises_off {SSM_NO_DEM_STO = 1 << 0, SSM_NO_WHITE_NOISE = 1 << 1, SSM_NO_DIFF = 1 << 2 }; //several noises can be turned off

enum ssm_print {SSM_PRINT_TRACE = 1 << 0, SSM_PRINT_X = 1 << 1, SSM_PRINT_HAT = 1 << 2, SSM_PRINT_PRED_RES = 1 << 3, SSM_PRINT_X_SMOOTH = 1 << 4, SSM_PRINT_ACC = 1 << 5, SSM_PIPE = 1 << 6, SSM_QUIET = 1 << 7, SSM_PRINT_COV = 1 << 8 };

typedef enum {SSM_SUCCESS = 1 << 0 , SSM_ERR_LIKE= 1 << 1, SSM_ERR_REM = 1 << 2} ssm_err_code;

#define SSM_BUFFER_SIZE (5000 * 1024)  /**< 5000 KB buffer size for settings.json inputs */
#define SSM_STR_BUFFSIZE 255 /**< buffer for log and error strings */
#define SSM_PATH_ROOT "./" /**< default root path for the results files (traces, ...) (has to be slash appended) */

#define SSM_WEB_APP 0 /**< webApp */

#define SSM_EPS_ABS 1e-6 /**< absolute error control for ODEs*/
#define SSM_EPS_REL 1e-3 /**< relative error control for ODEs*/

#define SSM_ZERO_LOG 1e-17 /**< smallest value that can be log transformed without being replaced by @c ZERO_LOG */
#define SSM_ONE_LOGIT 0.999999999 /**< largest value that can be logit transformed without being replaced by @c ONE_LOGIT */


/**
 * vector parameters in the user scale (all of them (infered or not))
 * ordered by: states, volatilities, process model parameters, observation model parameters
 */
typedef gsl_vector ssm_input_t;

/**
 * vector parameters in the natural scale (all of them (infered or not))
 * ordered as ssm_input_t
 */
typedef gsl_vector ssm_par_t;

/**
 * vector of *infered* parameters in the transformed scale (only the infered parameters)
 */
typedef gsl_vector ssm_theta_t;

/**
 * The variance-covariance matrix the transformed scale for the *infered* parameters (only the infered parameters)
 */
typedef gsl_matrix *ssm_var_t;




/**
 * the state variables (including including observed variables and
 * diffusions) and potientaly for kalman the covariance terms
 */
typedef struct  /* optionaly [N_DATA+1][J] for MIF and pMCMC "+1" is for initial condition (one time step before first data)  */
{
    gsl_vector *proj;    /**< values */

    double dt;           /**< the integration time step (for ODE solved with adaptive time step solvers) */
    double dt0;          /**< the integration time step initially picked by the user */
} ssm_X_t;


/**
 * The best estimates of the states variables (including observed
 * variables and diffusions) (weighted average of projected values,
 * weighted by the likelihood)
 */
typedef struct  /* ([N_DATA+1]) */
{
    int length;             /**< number of states */
    double *states;         /**< [self.length] best estimates */
    double **states_95;     /**< [self.length][2] 2.5% and 97.5% quantiles */
} ssm_X_hat_t;


/**
 * aggregation of parameters (e.g initial conditions, process model parameters...)
 */
typedef struct
{
    int length;           /**< number of parameters*/
    unsigned int *ind;    /**< [this.length] indexes of the parameters contained in the iterator */
} ssm_iterator_t;


/**
 * everything related to "parameters" (as opposed to states)
 */
typedef struct
{
    char *name; /**< name of the parameter */
    int order; /**< order of the parameter */

    double (*f) (double); /**< transformation (log, logit...) */
    double (*f_inv) (double); /**< inverse of f (f*f_inv=identity) */

    double (*f_derivative) (double); /**< derivative of f */
    double (*f_inv_derivative) (double); /**< derivative of f_inv */

    double (*prior) (double x); /**< prior */

    double (*f_par2user) (double); /**< from par to original user scale */
    double (*f_user2par) (double); /**< from original user scale to par */

} ssm_parameter_t;


/**
 * everything related to "states" including incidences, remainder...
 */
typedef struct
{
    char *name; /**< name of the state */
    int order; /**< order of the state */

    double (*f) (double); /**< transformation (log, logit...) */
    double (*f_inv) (double); /**< inverse of f (f*f_inv=identity) */

    double (*f_derivative) (double); /**< derivative of f */
    double (*f_inv_derivative) (double); /**< derivative of f_inv */

} ssm_state_t;


/**
 * navigating in the parameter / state space
 */
typedef struct
{
    //navigating withing par
    ssm_iterator_t *par_all;            /**< to iterate on every parameters */
    ssm_iterator_t *par_states;         /**< to iterate on the state variables *only* */
    ssm_iterator_t *par_remainders;     /**< to iterate on the remainders *only* */
    ssm_iterator_t *par_inc;            /**< to iterate on the state variables *only* */
    ssm_iterator_t *par_diff;           /**< to iterate on parameters following a diffusion *only* */
    ssm_iterator_t *par_noise;          /**< to iterate on white_noises sd *only* */

    //navigating within theta
    ssm_iterator_t *theta_all;                /**< to iterate on all the *infered* parameters */
    ssm_iterator_t *theta_no_states_no_diff;  /**< to iterate on the *infered* parameter of the process and observation models *not* following a diffusion */
    ssm_iterator_t *theta_states_diff;        /**< to iterate on the *infered* state variables *and* the *infered* parameters following a diffusion */

    int parameters_length; /**< total number of parameters (including non infered) but excluding covariate (present in ssm_calc_t) */
    ssm_parameter_t **parameters; /**< [this.parameters_length] <*/

    int states_length; /**< total number of states (including remainders and co) */
    ssm_state_t **states; /**< [this.states_length] <*/

} ssm_nav_t;


/**
 * likelihood values and associated quantities
 */
typedef struct
{
    int J;          /**< number of particles */
    int data_length;     /**< number of data points */
    int like_min;   /**< mimimun value of the likelihood */

    double least_square;

    double ess_n;               /**< effective sample size at n (sum(weight))^2 / sum(weight^2)*/
    double log_like_n ;         /**< log likelihood for the best parameter at n*/
    double log_like;            /**< log likelihood for the best parameter*/
    double *weights;            /**< [this.J] the weights */
    ssm_err_code *cum_status;  /**< [this.J] cumulated f_prediction status */

    unsigned int **select;      /**< [this.n_data][this.J] select is a vector with the indexes of the resampled particles. Note that we keep this.n_data values to keep genealogies */

    int n_all_fail;             /**< number of times when every particles had like < LIKE_MIN within one iteration */

    /* for bayesian methods */
    double log_like_prev;
    double log_like_new;

    /* prob priors */
    double *prior_probs;     /**< [this.J] prior probabilites */

} ssm_fitness_t;


/**
 * A row of data
 */
typedef struct { /* [n_data] */
    char *date;
    double *values;   /**< [n_ts] the data (values) one per time series */

    int ts_nonan_length; /**< number of time series without NaN at that time */
    unsigned int *ts_nonan; /**< [self.length] index of time series without NaN*/

    int states_reset_length; /**< number of states that must be reset to 0 */
    unsigned int *states_reset; /**< [self.length] index of states that must be reset to 0*/

} ssm_data_row_t;


/**
 * data
 */
typedef struct
{
    int n_ts;         /**< number of independent time series */
    int length;       /**< number of data points */

    int n_obs;       /**< the number of data point to taken into account for inference */

    char **names;     /**< [this.n_ts] name of the time series */

    ssm_data_row_t **rows; /**< [this.length] the data values */

    unsigned int *times;   /**< [this.length+1] ([0] + [times in days when the data were collected since the smallest t0]) */

    int length_nonan; /**< number of data points with at least one time series != NaN */
    unsigned int *ind_nonan; /**< [this.length_nonan] index of data points where there is at least one ts !=NaN */

} ssm_data_t;


/**
 * Everything needed to perform computations (possibly in parallel)
 * and store transiant states in a thread-safe way
 */
typedef struct /*[N_THREADS] : for parallel computing we need N_THREADS replication of the structure...*/
{
    int threads_length; /**< the total number of threads */
    int thread_id; /**< the id of the thread where the computation are being run */

    gsl_rng *randgsl; /**< random number generator */

    /////////////////
    //implementations
    /////////////////

    /* Euler multinomial */
    double **prob;      /**< [N_PAR_SV][number of output from the compartment]*/
    unsigned int **inc; /**< [N_PAR_SV][number of destinations] increments vector */

    /* Gillespie */
    //  double **reaction; /*reaction matrix*/

    /* ODE*/
    const gsl_odeiv2_step_type * T;
    gsl_odeiv2_control * control;
    gsl_odeiv2_step * step;
    gsl_odeiv2_evolve * evolve;
    gsl_odeiv2_system sys;
    double *yerr;

    /* SDE */
    double *y_pred; /**< used to store y predicted for Euler Maruyama */

    //multi-threaded sorting
    double *to_be_sorted;  /**< [J] array of the J particle to be sorted*/
    size_t *index_sorted;  /**< [J] index of the sorted weights used to compute 95% confidence interval */

    //interpolators for covariates
    gsl_interp_accel **acc;  /**< [N_PAR] an array of pointer to gsl_interp_accel */
    gsl_spline **spline;     /**< [N_PAR] an array of pointer to gsl_spline */

    /* references */
    ssm_theta_t *_par_natural; /**< Reference to the parameter is the
                                  natural scale (this.par_natural) used
                                  to pass it to some GSL function that
                                  only accept a *void. Such function
                                  only received ssm_calc_t * This
                                  reference should not be used outside
                                  from f_prediction_ functions. Outside
                                  these functions, this.par_natural is
                                  not guaranted to be defined. */

    ssm_nav_t *_par_natural; /**< ref to ssm_nav_t (same reason as ssm_theta_t) */
} ssm_calc_t;


/**
 * prediction function
 */
typedef ssm_err_code (*ssm_f_pred_t) (ssm_X_t *, double, double, ssm_theta_t *, ssm_nav_t *, struct ssm_calc_t *);


/**
 * observation function
 */
typedef ssm_err_code (*ssm_f_obs_t) (ssm_X_t *, double, double, ssm_theta_t *, ssm_nav_t *, struct ssm_calc_t *);


/**
 * Measuring durations
 */
typedef struct
{
    unsigned int d; ///< days
    unsigned int h; ///< hours
    unsigned int m; ///< minutes
    double s;       ///< seconds
} ssm_duration_t;


#endif
