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
#include <sysexits.h>

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

#include <jansson.h>

#include <zmq.h>
#include <pthread.h>

typedef enum {SSM_SMC = 1 << 0, SSM_MIF = 1 << 1, SSM_PMCMC = 1 << 2, SSM_KMCMC = 1 << 3, SSM_KALMAN = 1 << 4, SSM_KSIMPLEX = 1 << 5, SSM_SIMUL = 1 << 6, SSM_SIMPLEX = 1 << 7, SSM_WORKER = 1 << 8 } ssm_algo_t;
typedef enum {SSM_ODE, SSM_SDE, SSM_PSR, SSM_EKF} ssm_implementations_t;
typedef enum {SSM_NO_DEM_STO = 1 << 0, SSM_NO_WHITE_NOISE = 1 << 1, SSM_NO_DIFF = 1 << 2 } ssm_noises_off_t; //several noises can be turned off

typedef enum {SSM_PRINT_TRACE = 1 << 0, SSM_PRINT_X = 1 << 1, SSM_PRINT_HAT = 1 << 2, SSM_PRINT_DIAG = 1 << 3, SSM_PRINT_LOG = 1 << 4, SSM_PRINT_WARNING = 1 << 5 } ssm_print_t;

typedef enum {SSM_SUCCESS = 1 << 0 , SSM_ERR_LIKE= 1 << 1, SSM_ERR_REM_SV = 1 << 2, SSM_ERR_PRED = 1 << 3, SSM_ERR_KAL = 1 << 4, SSM_ERR_IC = 1 << 5, SSM_MH_REJECT = 1 << 6, SSM_ERR_PROPOSAL = 1 << 7, SSM_ERR_PRIOR = 1 << 8} ssm_err_code_t;

typedef enum {SSM_WORKER_J_PAR = 1 << 0, SSM_WORKER_D_X = 1 << 1, SSM_WORKER_FITNESS = 1 << 2 } ssm_worker_opt_t;

#define SSM_BUFFER_SIZE (10 * 1024)  /**< 1000 KB buffer size */
#define SSM_STR_BUFFSIZE 255 /**< buffer for log and error strings */


#define SSM_WEB_APP 0 /**< webApp */

#define SSM_ZERO_LOG 1e-17 /**< smallest value that can be log transformed without being replaced by @c ZERO_LOG */
#define SSM_ONE_LOGIT 0.999999999 /**< largest value that can be logit transformed without being replaced by @c ONE_LOGIT */


/**
 * vector parameters in the user scale (all of them (infered or not))
 * ordered by: states, volatilities, noises, process model parameters, observation model parameters
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
typedef gsl_matrix ssm_var_t;


typedef struct _nav ssm_nav_t;

/**
 * Everything needed to perform computations (possibly in parallel)
 * and store transiant states in a thread-safe way
 */
typedef struct ssm_calc_t /*[N_THREADS] : for parallel computing we need N_THREADS replication of the structure...*/
{
    int seed;
    int threads_length; /**< the total number of threads */
    int thread_id;      /**< the id of the thread where the computation are being run */

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
    const gsl_odeiv2_step_type *T;
    gsl_odeiv2_control *control;
    gsl_odeiv2_step *step;
    gsl_odeiv2_evolve *evolve;
    gsl_odeiv2_system sys;
    double *yerr;

    /* SDE */
    double *y_pred; /**< used to store y predicted for Euler Maruyama */

    /* Kalman */
    void (*eval_Q)(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, struct ssm_calc_t *calc);

    gsl_vector *_pred_error;    /**< [nav->observed_length] */
    gsl_vector *_zero;          /**< [nav->observed_length] */
    gsl_matrix *_St;            /**< [nav->observed_length][nav->observed_length] */
    gsl_matrix *_Stm1;          /**< [nav->observed_length][nav->observed_length] */
    gsl_matrix *_Rt;            /**< [nav->observed_length][nav->observed_length] */
    gsl_matrix *_Ht;            /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->observed_length] */
    gsl_matrix *_Kt;            /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->observed_length] */
    gsl_matrix *_Tmp_N_SV_N_TS; /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->observed_length] */
    gsl_matrix *_Tmp_N_TS_N_SV; /**< [nav->observed_length][nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_matrix *_Q;             /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_matrix *_FtCt;          /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_matrix *_Ft;            /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_vector *_eval;      /**< [nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_matrix *_evec;      /**< [nav->states_sv_inc->length + nav->states_diff->length][nav->states_sv_inc->length + nav->states_diff->length] */
    gsl_eigen_symmv_workspace *_w_eigen_vv;  /**< workspace to compute eigen values and eigen vector for symmetric matrix */

    //multi-threaded sorting
    int J;                 /**< ssm_fitness_t->J */
    double *to_be_sorted;  /**< [this->J] array of the J particle to be sorted*/
    size_t *index_sorted;  /**< [this->J] index of the sorted weights used to compute 95% confidence interval */

    //interpolators for covariates
    int covariates_length;   /**< number of covariates */
    gsl_interp_accel **acc;  /**< [self.covariates_length] an array of pointer to gsl_interp_accel */
    gsl_spline **spline;     /**< [self.covariates_length] an array of pointer to gsl_spline */

    /* references */
    ssm_par_t *_par; /**< Reference to the parameter is the natural
                        scale (this.par) used to pass it to
                        some GSL function that only accept a
                        *void. Such function only received ssm_calc_t
                        * This reference should not be used outside
                        from f_prediction_ functions. Outside these
                        functions, this.par_natural is not guaranted
                        to be defined. */

    ssm_nav_t *_nav; /**< ref to ssm_nav_t (same reason as ssm_par_t) */
} ssm_calc_t;



/**
 * the state variables (including including observed variables and
 * diffusions) and potientaly for kalman the covariance terms
 */
typedef struct  /* optionaly [N_DATA+1][J] for MIF and pMCMC "+1" is for initial condition (one time step before first data)  */
{
    int length;
    double *proj; /**< [this.length] values */

    double dt;  /**< the integration time step (for ODE solved with adaptive time step solvers) */
    double dt0; /**< the integration time step initially picked by the user */

} ssm_X_t;


/**
 * The best estimates of the states variables (including observed
 * variables and diffusions and remainders) (weighted average of projected values,
 * weighted by the likelihood)
 */
typedef struct  /* ([N_DATA+1]) */
{
    int states_length;      /**< number of states (== nav->states_sv_inc->length + nav->states_diff->length) */
    double *states;         /**< [this.length] best estimates */
    double **states_95;     /**< [this.length][2] 2.5% and 97.5% quantiles */

    int remainders_length;
    double *remainders;     /**< [this.remainders_length] best estimates */
    double **remainders_95; /**< [this.remainders_length][2] 2.5% and 97.5% quantiles */

    int observed_length;    /**< number of observed variable (== ssm_nav_t->observed_length) */
    double *observed;       /**< [this.observed_length] best estimates */
    double **observed_95;   /**< [this.observed_length][2] 2.5% and 97.5% quantiles */

} ssm_hat_t;



/**
 * everything related to "parameters" (as opposed to states)
 */
typedef struct
{
    char *name; /**< name of the parameter */
    int offset; /**< order of the parameter in ssm_par_t and ssm_input_t */
    int offset_theta; /**< order of the parameter in ssm_theta_t (-1 if the parameter is not in ssm_theta_t)*/

    double (*f) (double x_input); /**< transformation (log, logit...) */
    double (*f_inv) (double x_theta); /**< inverse of f (f*f_inv=identity) */

    double (*f_der) (double x_theta); /**< derivative of f */
    double (*f_der_inv) (double x_theta); /**< derivative of f_inv */
    double (*f_der2_inv) (double x_theta); /**< 2nd derivative of f_inv */

    double (*f_prior) (double x_input); /**< prior */

    double (*f_user2par) (double x_input, ssm_input_t *, ssm_calc_t *); /**< from original user scale to par */

    double (*f_2prior) (double x_input, ssm_hat_t *hat, ssm_par_t *par, ssm_calc_t *calc, double t); /**< Transform X coming from hat into in input in the user unit (prior (resource)). This is usefull when we want to pipe the end of a simulation. In this case, x_input will be ignored and reconstructed from hat and par using the to_resource function provided by the user. x_input is provided for parameters different from states initial condition. In this case x_input is simply returned */

} ssm_parameter_t;



/**
 * aggregation of parameters (e.g initial conditions, process model parameters...)
 */
typedef struct
{
    int length;           /**< number of parameters*/
    ssm_parameter_t **p;  /**< [this.length] */
} ssm_it_parameters_t;




/**
 * everything related to "states" including incidences, remainder...
 */
typedef struct
{
    char *name; /**< name of the state */
    int offset; /**< order of the state in ssm_X_t or ssm_hat_t (for remainders it is hat->remainders) */

    ssm_parameter_t *ic;  /**< pointer to the initial condition (or NULL) */

    double (*f) (double x_X);     /**< transformation (log, logit...) */
    double (*f_inv) (double x_X); /**< inverse of f (f*f_inv=identity) */

    double (*f_der) (double x_X);     /**< derivative of f */
    double (*f_der_inv) (double x_X); /**< derivative of f_inv */
    double (*f_der2_inv) (double x_X); /**< 2nd derivative of f_inv */

    double (*f_remainder) (ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, double t); /**< compute the remainder value */
    double (*f_remainder_var) (ssm_X_t *X, ssm_calc_t *calc, ssm_nav_t *nav, double t); /**< compute the remainder value */
} ssm_state_t;


/**
 * aggregation of states (state variables, remainder...)
 */
typedef struct
{
    int length;           /**< number of parameters*/
    ssm_state_t **p;  /**< [this.length] */
} ssm_it_states_t;


/**
 * everything related to observed variable
 */
typedef struct
{
    char *name; /**< name of the observed variable */
    int offset; /**< order of the observed variable in nav.observed */

    double (*f_likelihood) (double y, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, double t);
    double (*f_obs_mean)             (ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, double t);
    double (*f_obs_var)              (ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, double t);
    double (*f_obs_ran)              (ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, double t);
    double (*f_var_pred)             (ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, double t);
} ssm_observed_t;


/**
 * navigating in the parameter / state space
 */
struct _nav
{
    ssm_implementations_t implementation;
    ssm_noises_off_t noises_off;
    ssm_print_t print;


    FILE *X;
    FILE *hat;
    FILE *diag;
    FILE *trace;

    int parameters_length;         /**< total number of parameters (including non infered) but excluding covariate (present in ssm_calc_t) */
    ssm_parameter_t **parameters; /**< [this.parameters_length] <*/

    int states_length;    /**< total number of states (including remainders and co) */
    ssm_state_t **states; /**< [this.states_length] <*/

    int observed_length;       /**< total number of observed variables */
    ssm_observed_t **observed; /**< [this.observed_length] */

    //navigating withing par
    ssm_it_states_t *states_sv;         /**< to iterate on the state variables (not including remainders or incidences) *only* */
    ssm_it_states_t *states_remainders; /**< to iterate on the remainders *only* */
    ssm_it_states_t *states_inc;        /**< to iterate on the incidences *only* */
    ssm_it_states_t *states_sv_inc;     /**< to iterate on the state variables and incidences *only* */
    ssm_it_states_t *states_diff;       /**< to iterate on states following a diffusion *only* */

    ssm_it_parameters_t *par_all;       /**< to iterate on every parameters */
    ssm_it_parameters_t *par_noise;     /**< to iterate on white_noises sd *only* */
    ssm_it_parameters_t *par_disp;      /**< to iterate on parameter only involve in dispersion *only* */
    ssm_it_parameters_t *par_icsv;      /**< to iterate on the initial condition of the state variables *only* */
    ssm_it_parameters_t *par_icdiff;    /**< to iterate on the initial condition of the diffusions *only* */

    //navigating within theta
    ssm_it_parameters_t *theta_all;                /**< to iterate on all the *infered* parameters */
    ssm_it_parameters_t *theta_no_icsv_no_icdiff;  /**< to iterate on the *infered* parameter of the process and observation models *not* being initial conditions of state variables or initial conditions of diffusions */
    ssm_it_parameters_t *theta_icsv_icdiff;        /**< to iterate on the *infered* initial conditions of the state variables *and* the *infered* initial conditions of the parameters following a diffusion */
};


/**
 * likelihood values and associated quantities
 */
typedef struct
{
    int J;           /**< number of particles */
    int data_length; /**< number of data points */
    double like_min; /**< mimimun value of the likelihood */
    double log_like_min; /**< mimimun value of the log likelihood */

    double ess_n;               /**< effective sample size at n (sum(weight))^2 / sum(weight^2)*/
    double log_like_n ;         /**< log likelihood for the best parameter at n*/
    double log_like;            /**< log likelihood for the best parameter*/

    double *weights;            /**< [this.J] the weights */
    unsigned int **select;      /**< [this.data_length][this.J] select is a vector with the indexes of the resampled particles. Note that we keep this.n_data values to keep genealogies */

    ssm_err_code_t *cum_status;   /**< [this.J] cumulated f_prediction status */

    int n_all_fail;             /**< number of times when every particles had like < LIKE_MIN within one iteration */

    /* for bayesian methods */
    double log_like_prev;

    double log_prior;
    double log_prior_prev;

    //summary quantities
    double AIC;
    double AICc;
    double DIC;
    double summary_log_likelihood; /**< ALWAYS without prior! the final log likelihood for non MCMC methods and the best likelihood found during the trace for MCMC methods*/
    double summary_log_ltp; /**< log of likelifhood times prior (same as above, for MCMC methods this is the best value found during the trace) */
    double summary_sum_squares;
    int n; /**< number of data point non NA consumed */

    double _min_deviance; /**< min(deviance)*/
    double _deviance_cum; /**< sum(deviance)*/

} ssm_fitness_t;


/**
 * A row of data
 */
typedef struct { /* [n_data] */
    char *date;
    unsigned int time;          /**< times in days since the smallest dates_t0 when the data were collected */

    int ts_nonan_length;        /**< number of time series without NaN at that time */
    ssm_observed_t **observed;  /**< [self.ts_nonan_length] */
    double *values;             /**< [ts_nonan_length] the non NAN data (values) one per time series */

    int states_reset_length;    /**< number of states that must be reset to 0 */
    ssm_state_t **states_reset; /**< [self.states_reset_length] array of pointer to the states to reset */

} ssm_row_t;


/**
 * data
 */
typedef struct
{
    int length;              /**< number of rows (aligned data points) */
    int ts_length;           /**< the number of time series */
    char* date_t0;           /**< the smallest date (ISO 8601) before the first data point*/
    int n_obs;               /**< the number of rows (aligned data points) to taken into account for inference */
    int n_obs_nonan;         /**< the number of rows (aligned data points) to taken into account for inference discarding lines where all ts are NaN  */

    ssm_row_t **rows;   /**< [this.length] the data values */
    int length_nonan;        /**< number of data points with at least one time series != NaN */
    unsigned int *ind_nonan; /**< [this.length_nonan] index of data points where there is at least one ts !=NaN */

} ssm_data_t;



/**
 * prediction function
 */
typedef ssm_err_code_t (*ssm_f_pred_t) (ssm_X_t *, double, double, ssm_par_t *, ssm_nav_t *, ssm_calc_t *);


/**
 * options
 */
typedef struct
{
    ssm_algo_t algo;
    ssm_algo_t worker_algo;

    ssm_implementations_t implementation;
    ssm_noises_off_t noises_off;
    ssm_print_t print;

    int id;                  /**< unique integer identifier that will be used as seed and be appended to the output files */
    int flag_seed_time;      /**< seed with the local time ((unsigned) time(NULL)) */
    int flag_prior;          /**< add log(prior) to the estimated log likelihood */
    double dt;               /**< integration time step in days */
    double eps_abs;          /**< absolute error for adaptive step-size control */
    double eps_rel;          /**< relative error for adaptive step-size control */
    char *freeze_forcing;    /**< freeze the metadata to their value at the specified ISO8601 date */
    char *root;              /**< root path where the outputs will be stored */
    char *next;              /**< write the outputed parameters in a file prefixed by the argument */
    int n_thread;            /**< number of threads */
    double like_min;         /**< particles with likelihood smaller that like_min are considered lost */
    int J;                   /**< number of particles */
    int n_obs;               /**< number of observations to be fitted (for tempering) */
    char *interpolator;      /**< gsl interpolator for metadata */
    int n_iter;              /**< number of iterations */
    double a;                /**< cooling factor (scales standard deviation) */
    double b;                /**< re-heating (inflation) (scales standard deviation of the proposal) */
    double L;                /**< lag for fixed lag smoothing (proportion of the data) */
    int m_switch;            /**< iteration number when we switch (for initial covariance to empirical or from different update formulae) */
    int flag_ic_only;        /**< only fit the initial condition using fixed lag smoothing */
    int eps_switch;          /**< select number of burnin iterations before tuning epsilon */
    double eps_max;          /**< maximum value allowed for epsilon */
    double flag_smooth;      /**< tune epsilon with the value of the acceptance rate obtained with exponential smoothing */
    double alpha;            /**< smoothing factor of exponential smoothing used to compute the smoothed acceptance rate (low values increase degree of smoothing) */
    int n_traj;              /**< number of trajectories stored */
    int flag_tcp;            /**< dispatch particles across machines */
    int flag_least_squares;  /**< optimize the sum of square instead of the likelihood */
    double size_stop;        /**< simplex size used as a stopping criteria */
    int freq;                /**< print the outputs (and reset incidences to 0 if any) every specified days */
    char *start;             /**< ISO 8601 date when the simulation starts*/
    char *end;               /**< ISO 8601 date when the simulation ends*/
    char *server;            /**< domain name or IP address of the particule server (e.g 127.0.0.1) */
    int flag_no_filter;      /**< do not filter */
} ssm_options_t;




/**
 * Adaptive tunning of MCMC algo
 */
typedef struct
{
    double ar;          /**< the global acceptance rate */
    double ar_smoothed; /**< as computed with exponential smoothing (http://en.wikipedia.org/wiki/Exponential_smoothing) */

    double eps;         /**< epsilon factor */
    double eps_a;       /**< cooling factor */
    double eps_max;     /**< max value for epsilon */
    int    eps_switch;  /**< number of iterations before tuning epsilon */


    int m_switch;       /**< number of iterations using empirical covariance */

    int flag_smooth; /**<boolean: do we tune epsilon with the value of the acceptance rate obtained with exponential smoothing ? (1 yes 0 no) */
    double alpha; /**< smoothing factor (the term smoothing factor is
                     something of a misnomer, as larger values of
                     alpha actually reduce the level of smoothing, and
                     in the limiting case with alpha = 1 the output
                     series is just the same as the original series
                     (with lag of one time unit). */

    double *mean_sampling;         /**< [ssm_nav_t->theta_all->length] Em(X) 1st order mean needed to compute the sampling covariance */
    gsl_matrix *var_sampling;      /**< [ssm_nav_t->theta_all->length][ssm_nav_t->theta_all->length] Sampling covariance */

} ssm_adapt_t;



typedef struct
{
    int id;
    void *context; /**< zmq context */
    ssm_worker_opt_t wopts;
    int J_chunk;
    ssm_data_t *data;
    ssm_par_t **J_par;
    ssm_X_t ***D_J_X;
    ssm_calc_t **calc;
    ssm_nav_t *nav;
    ssm_fitness_t *fitness;
    ssm_f_pred_t f_pred;
} ssm_params_worker_inproc_t;


typedef struct 
{
    int flag_tcp;
    int inproc_length; /**< number of inproc worker */
    ssm_worker_opt_t wopts;

    void *context;
    void *sender;
    void *receiver;
    void *controller;
    ssm_params_worker_inproc_t *params;

    pthread_t *workers;
} ssm_workers_t;


/****************************/
/* core function signatures */
/****************************/

/* alloc_c.c */
char *ssm_c1_new(int n);
char **ssm_c2_new(int n, int p);
void ssm_c2_free(char **tab, int n);

/* alloc_d.c */
double *ssm_d1_new(int n);
double **ssm_d2_new(int n, int p);
void ssm_d2_free(double **tab, int n);
double ***ssm_d3_new(int n, int p1, int p2);
void ssm_d3_free(double ***tab, int n, int p1);
double ****ssm_d4_new(int n, int p1, int p2, int p3);
void ssm_d4_free(double ****tab, int n, int p1, int p2);
double **ssm_d2_var_new(int n, unsigned int *p);
double ***ssm_d3_var_new(int n, unsigned int *p1, unsigned int **p2);
void ssm_d3_var_free(double ***tab, int n, unsigned int *p1);
double ***ssm_d3_varp1_new(int n, unsigned int *p1, int p2);
double ***ssm_d3_varp2_new(int n, unsigned int p1, unsigned int *p2);

/* alloc_i.c */
int *ssm_i1_new(int n);

/* alloc_st.c */
size_t *ssm_st1_new(int n);
size_t **ssm_st2_new(int n, int p);
void ssm_st2_free(size_t **tab, int n);

/* alloc_u.c */
unsigned int *ssm_u1_new(int n);
unsigned int **ssm_u2_new(int n, int p);
void ssm_u2_free(unsigned int **tab, int n);
unsigned int ***ssm_u3_new(int n, int p1, int p2);
void ssm_u3_free(unsigned int ***tab, int n, int p1);
unsigned int ****ssm_u4_new(int n, int p1, int p2, int p3);
void ssm_u4_free(unsigned int ****tab, int n, int p1, int p2);
unsigned int **ssm_u2_var_new(int n, unsigned int *p);
unsigned int ***ssm_u3_var_new(int n, unsigned int *p1, unsigned int **p2);
void ssm_u3_var_free(unsigned int ***tab, int n, unsigned int *p1);
unsigned int ***ssm_u3_varp1_new(int n, unsigned int *p1, int p2);
unsigned int ***ssm_u3_varp2_new(int n, unsigned int p1, unsigned int *p2);

/* build.c */
void ssm_input_free(ssm_input_t *input);
ssm_par_t *ssm_par_new(ssm_input_t *input, ssm_calc_t *calc, ssm_nav_t *nav);
void ssm_par_free(ssm_par_t *par);
ssm_theta_t *ssm_theta_new(ssm_input_t *input, ssm_nav_t *nav);
void ssm_theta_free(ssm_theta_t *theta);
ssm_var_t *ssm_var_new(json_t *jparameters, ssm_nav_t *nav);
void ssm_var_free(ssm_var_t *var);
ssm_it_states_t *_ssm_it_states_new(int length);
void _ssm_it_states_free(ssm_it_states_t *it);
ssm_it_parameters_t *_ssm_it_parameters_new(int length);
void _ssm_it_parameters_free(ssm_it_parameters_t *it);
ssm_nav_t *ssm_nav_new(json_t *jparameters, ssm_options_t *opts);
void _ssm_observed_free(ssm_observed_t *observed);
void _ssm_parameter_free(ssm_parameter_t *parameter);
void _ssm_state_free(ssm_state_t *state);
void ssm_nav_free(ssm_nav_t *nav);
ssm_data_t *ssm_data_new(json_t *jdata, ssm_nav_t *nav, ssm_options_t *opts);
void _ssm_row_free(ssm_row_t *row);
void ssm_data_free(ssm_data_t *data);
void ssm_data_adapt_to_simul(ssm_data_t *data, json_t *jdata, ssm_nav_t *nav, ssm_options_t *opts);
ssm_calc_t *ssm_calc_new(json_t *jdata, ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, ssm_options_t *opts, int thread_id);
void ssm_calc_free(ssm_calc_t *calc, ssm_nav_t *nav);
ssm_calc_t **ssm_N_calc_new(json_t *jdata, ssm_nav_t *nav, ssm_data_t *data, ssm_fitness_t *fitness, ssm_options_t *opts);
void ssm_N_calc_free(ssm_calc_t **calc, ssm_nav_t *nav);
ssm_options_t *ssm_options_new(void);
void ssm_options_free(ssm_options_t *opts);
ssm_fitness_t *ssm_fitness_new(ssm_data_t *data, ssm_options_t *opts);
void ssm_fitness_free(ssm_fitness_t *fitness);
int _ssm_dim_X(ssm_nav_t *nav);
ssm_X_t *ssm_X_new(ssm_nav_t *nav, ssm_options_t *opts);
void ssm_X_free(ssm_X_t *X);
ssm_X_t **ssm_J_X_new(ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_options_t *opts);
void ssm_J_X_free(ssm_X_t **X, ssm_fitness_t *fitness);
ssm_X_t **ssm_D_X_new(ssm_data_t *data, ssm_nav_t *nav, ssm_options_t *opts);
void ssm_D_X_free(ssm_X_t **X, ssm_data_t *data);
ssm_X_t ***ssm_D_J_X_new(ssm_data_t *data, ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_options_t *opts);
void ssm_D_J_X_free(ssm_X_t ***X, ssm_data_t *data, ssm_fitness_t *fitness);
ssm_hat_t *ssm_hat_new(ssm_nav_t *nav);
void ssm_hat_free(ssm_hat_t *hat);
ssm_hat_t **ssm_D_hat_new(ssm_data_t *data, ssm_nav_t *nav);
void ssm_D_hat_free(ssm_hat_t **hat, ssm_data_t *data);
ssm_adapt_t *ssm_adapt_new(ssm_nav_t *nav, ssm_options_t * opts);
void ssm_adapt_free(ssm_adapt_t *adapt);

/* load.c */
json_t *ssm_load_json_stream(FILE *stream);
json_t *ssm_load_json_file(const char *path);
json_t *ssm_load_data(ssm_options_t *opts);
void ssm_theta2input(ssm_input_t *input, ssm_theta_t *theta, ssm_nav_t *nav);
void ssm_input2par(ssm_par_t *par, ssm_input_t *input, ssm_calc_t *calc, ssm_nav_t *nav);
void ssm_par2X(ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav);
void ssm_mcmc_results2X(ssm_X_t *X, json_t *jprediction_j, ssm_calc_t *calc, ssm_nav_t *nav);
unsigned int *ssm_load_ju1_new(json_t *container, char *name);
double *ssm_load_jd1_new(json_t *container, char *name);
char **ssm_load_jc1_new(json_t *container, const char *name);

/* options.c */
void ssm_options_load(ssm_options_t *opts, ssm_algo_t algo, int argc, char *argv[]);
void ssm_options_set_implementation(ssm_options_t *opts, ssm_algo_t algo, int argc, char *argv[]);

/* fitness.c */
double ssm_sanitize_log_likelihood(double log_like, ssm_row_t *row, ssm_fitness_t *fitness, ssm_nav_t *nav);
double ssm_log_likelihood(ssm_row_t *row, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *fitness);
double ssm_sum_square(ssm_row_t *row, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_fitness_t *fitness);
void ssm_aic(ssm_fitness_t *fitness, ssm_nav_t *nav, double log_like);
void ssm_dic_init(ssm_fitness_t *fitness, double log_like, double log_prior);
void ssm_dic_update(ssm_fitness_t *fitness, double log_like, double log_prior);
void ssm_dic_end(ssm_fitness_t *fitness, ssm_nav_t *nav, int m);

/* mvn.c */
int ssm_rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, double sd_fac, gsl_vector *result);
double ssm_dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, double sd_fac);
void ssm_adapt_var(ssm_adapt_t *adapt, ssm_theta_t *x, int m);

/* prediction_util.c */
void ssm_X_copy(ssm_X_t *dest, ssm_X_t *src);
void ssm_X_reset_inc(ssm_X_t *X, ssm_row_t *row, ssm_nav_t *nav);
void ssm_ran_multinomial (const gsl_rng * r, const size_t K, unsigned int N, const double p[], unsigned int n[]);
double ssm_correct_rate(double rate, double dt);
ssm_err_code_t ssm_check_no_neg_sv_or_remainder(ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, double t);
ssm_f_pred_t ssm_get_f_pred(ssm_nav_t *nav);
ssm_err_code_t ssm_f_prediction_ode                           (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_dem_sto_no_white_noise (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_dem_sto_no_diff        (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_white_noise_no_diff    (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_dem_sto                (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_white_noise            (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_no_diff                   (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_sde_full                      (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_psr                           (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
ssm_err_code_t ssm_f_prediction_psr_no_diff                   (ssm_X_t *p_X, double t0, double t1, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* smc.c */
int ssm_weight(ssm_fitness_t *fitness, ssm_row_t *row, ssm_nav_t *nav, int n);
void ssm_systematic_sampling(ssm_fitness_t *fitness, ssm_calc_t *calc, int n);
void ssm_resample_X(ssm_fitness_t *fitness, ssm_X_t ***J_p_X, ssm_X_t ***J_p_X_tmp, int n);
void ssm_swap_X(ssm_X_t ***X, ssm_X_t ***tmp_X);

/* transform.c */
double ssm_f_id(double x);
double ssm_f_der_id(double x);
double ssm_f_der2_inv_id(double x);
double ssm_f_log(double x);
double ssm_f_inv_log(double x);
double ssm_f_logit(double x);
double ssm_f_inv_logit(double x);
double ssm_f_logit_ab(double x, double a, double b);
double ssm_f_inv_logit_ab(double x, double a, double b);
double ssm_f_der_log(double x);
double ssm_f_der_inv_log(double x);
double ssm_f_der2_inv_log(double x);
double ssm_f_der_logit(double x);
double ssm_f_der_inv_logit(double x);
double ssm_f_der2_inv_logit(double x);
double ssm_f_der_logit_ab(double x, double a, double b);
double ssm_f_der_inv_logit_ab(double x, double a, double b);
double ssm_f_der2_inv_logit_ab(double x, double a, double b);
double ssm_f_user_par_id(double x, ssm_input_t *par, ssm_calc_t *calc);
double ssm_f_2prior_id(double x, ssm_hat_t *hat, ssm_par_t *par, ssm_calc_t *calc, double t);

/* util.c */
int ssm_in_par(ssm_it_parameters_t *it, const char *name);
int ssm_in_jarray(json_t *array, const char *name);
const gsl_interp_type *ssm_str_to_interp_type(const char *optarg);
int ssm_sanitize_n_threads(int n_threads, ssm_fitness_t *fitness);

/* print.c */
void ssm_print_log(char *msg);
void ssm_print_warning(char *msg);
void ssm_print_err(char *msg);
void ssm_json_dumpf(FILE *stream, const char *flag, json_t *msg);
void ssm_pipe_theta(FILE *stream, json_t *jparameters, ssm_theta_t *theta, ssm_var_t *var, ssm_fitness_t *fitness, ssm_nav_t *nav, ssm_options_t *opts);
void ssm_pipe_hat(FILE *stream, json_t *jparameters, ssm_input_t *input, ssm_hat_t *hat, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav, ssm_options_t *opts, double t);
void ssm_print_header_X(FILE *stream, ssm_nav_t *nav);
void ssm_print_X(FILE *stream, ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_row_t *row, const int index);
void ssm_print_header_trace(FILE *stream, ssm_nav_t *nav);
void ssm_print_trace(FILE *stream, ssm_theta_t *theta, ssm_nav_t *nav, const double fitness, const int index);
void ssm_print_header_pred_res(FILE *stream, ssm_nav_t *nav);
void ssm_print_pred_res(FILE *stream, ssm_X_t **J_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_data_t *data, ssm_row_t *row, ssm_fitness_t *fitness);
void ssm_print_header_hat(FILE *stream, ssm_nav_t *nav);
void ssm_print_hat(FILE *stream, ssm_hat_t *hat, ssm_nav_t *nav, ssm_row_t *row);
void ssm_sample_traj_print(FILE *stream, ssm_X_t ***D_J_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness, const int index);
int ssm_sample_traj_print2(FILE *stream, ssm_X_t ***D_J_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness, const int index);
void ssm_print_header_ar(FILE *stream);
void ssm_print_ar(FILE *stream, ssm_adapt_t *adapt, const int index);

/* hat.c */
void ssm_ci95(double *hat_95, ssm_calc_t *calc, ssm_fitness_t *fitness);
void ssm_hat_eval(ssm_hat_t *hat, ssm_X_t **J_X, ssm_par_t **J_par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_fitness_t *fitness, const double t, int is_J_par);


/* bayes.c */
ssm_err_code_t ssm_log_prob_proposal(double *log_proposal, ssm_theta_t *proposed, ssm_theta_t *theta, ssm_var_t *var, double sd_fac, ssm_nav_t *nav, int is_mvn);
ssm_err_code_t ssm_log_prob_prior(double *log_prior, ssm_theta_t *theta, ssm_nav_t *nav, ssm_fitness_t *fitness);
ssm_err_code_t ssm_metropolis_hastings(ssm_fitness_t *fitness, double *alpha, ssm_theta_t *proposed, ssm_theta_t *theta, gsl_matrix *var, double sd_fac , ssm_nav_t *nav, ssm_calc_t *calc, int is_mvn);
ssm_var_t *ssm_adapt_eps_var_sd_fac(double *sd_fac, ssm_adapt_t *a, ssm_var_t *var, ssm_nav_t *nav, int m);
void ssm_adapt_ar(ssm_adapt_t *a, int is_accepted, int m);
void ssm_theta_ran(ssm_theta_t *proposed, ssm_theta_t *theta, ssm_var_t *var, double sd_fac, ssm_calc_t *calc, ssm_nav_t *nav, int is_mvn);
int ssm_theta_copy(ssm_theta_t *dest, ssm_theta_t *src);
int ssm_par_copy(ssm_par_t *dest, ssm_par_t *src);
void ssm_sample_traj(ssm_X_t **D_X, ssm_X_t ***D_J_X, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness);
void ssm_sample_traj2(ssm_X_t **D_X, ssm_X_t ***D_J_X, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness, const int j_select);

/* simplex.c */
double ssm_simplex(ssm_theta_t *theta, ssm_var_t *var, void *params, double (*f_simplex)(const gsl_vector *x, void *params), ssm_nav_t *nav, ssm_options_t *opts);

/* workers.c */
void *ssm_worker_inproc(void *params);
ssm_workers_t *ssm_workers_start(ssm_X_t ***D_J_X, ssm_par_t **J_par, ssm_data_t *data, ssm_calc_t **calc, ssm_fitness_t *fitness, ssm_f_pred_t f_pred, ssm_nav_t *nav, ssm_options_t *opts, ssm_worker_opt_t wopts);
void ssm_workers_stop(ssm_workers_t *workers);

/* special functions */
double heaviside(double x);
double ramp(double x);
double slowstep(double x, double d);
double sigmoid(double x, double shape, double shift);
 
/******************************/
/* kalman function signatures */
/******************************/

/* kalman/ekf.c */
ssm_err_code_t _ssm_check_and_correct_Ct(ssm_X_t *X, ssm_calc_t *calc, ssm_nav_t *nav);
ssm_err_code_t ssm_kalman_gain_computation(ssm_row_t *row, double t, ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav);
ssm_err_code_t ssm_kalman_update(ssm_fitness_t *fitness, ssm_X_t *X, ssm_row_t *row, double t, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav);
double ssm_diff_derivative(double jac_tpl, const double X[], ssm_state_t *state);
void ssm_kalman_reset_Ct(ssm_X_t *X, ssm_nav_t *nav);

/******************************/
/* mif function signatures */
/******************************/

/* mif/mif_util.c */
double ssm_mif_cooling(ssm_options_t *opts, int m);
void ssm_mif_scale_var(ssm_var_t *var, ssm_data_t *data, ssm_nav_t *nav);
void ssm_mif_patch_like_prior(double *like, ssm_fitness_t *fitness, ssm_theta_t **J_theta, ssm_data_t *data, ssm_nav_t *nav, const int n, const int lag);
void ssm_mif_mean_var_theta_theoretical(double *theta_bart, double *theta_Vt, ssm_theta_t **J_theta, ssm_var_t *var, ssm_fitness_t *fitness, ssm_nav_t *nav, double var_fac);
void ssm_mif_resample_and_mutate_theta(ssm_fitness_t *fitness, ssm_theta_t **J_theta, ssm_theta_t **J_theta_tmp, ssm_var_t *var, ssm_calc_t **calc, ssm_nav_t *nav, double sd_fac, int n);
void ssm_mif_fixed_lag_smoothing(ssm_theta_t *mle, ssm_theta_t **J_theta, ssm_fitness_t *fitness, ssm_nav_t *nav);
void ssm_mif_update_average(ssm_theta_t *mle, double **D_theta_bart, ssm_data_t *data, ssm_nav_t *nav);
void ssm_mif_update_ionides(ssm_theta_t *mle, ssm_var_t *var, double **D_theta_bart, double **D_theta_Vt, ssm_data_t *data, ssm_nav_t *nav, ssm_options_t *opts, double cooling);
void ssm_mif_print_header_mean_var_theoretical_ess(FILE *stream, ssm_nav_t *nav);
void ssm_mif_print_mean_var_theoretical_ess(FILE *stream, double *theta_bart, double *theta_Vt, ssm_fitness_t *fitness, ssm_nav_t *nav , ssm_row_t *row, int m);

/******************************/
/* worker function signatures */
/******************************/

/* worker/worker_util.c */
void ssm_zmq_send_par(void *socket, ssm_par_t *par, int zmq_options);
void ssm_zmq_recv_par(ssm_par_t *par, void *socket);
void ssm_zmq_send_X(void *socket, ssm_X_t *X, int zmq_options);
void ssm_zmq_recv_X(ssm_X_t *X, void *socket);

/*********************************/
/* templated function signatures */
/*********************************/

/* input_template.c */
ssm_input_t *ssm_input_new(json_t *jparameters, ssm_nav_t *nav);

/* transform_template.c */
ssm_parameter_t **_ssm_parameters_new(int *parameters_length);
ssm_state_t **_ssm_states_new(int *states_length, ssm_parameter_t **parameters);

/* check_IC_template */
ssm_err_code_t ssm_check_ic(ssm_par_t *par, ssm_calc_t *calc);

/* iterator_template.c */
ssm_it_states_t *ssm_it_states_sv_new(ssm_state_t **states);
ssm_it_states_t *ssm_it_states_remainders_new(ssm_state_t **states);
ssm_it_states_t *ssm_it_states_inc_new(ssm_state_t **states);
ssm_it_states_t *ssm_it_states_sv_inc_new(ssm_state_t **states);
ssm_it_states_t *ssm_it_states_diff_new(ssm_state_t **states);
ssm_it_parameters_t *ssm_it_parameters_all_new(ssm_parameter_t **parameters);
ssm_it_parameters_t *ssm_it_parameters_noise_new(ssm_parameter_t **parameters);
ssm_it_parameters_t *ssm_it_parameters_disp_new(ssm_parameter_t **parameters);
ssm_it_parameters_t *ssm_it_parameters_vol_new(ssm_parameter_t **parameters);
ssm_it_parameters_t *ssm_it_parameters_icsv_new(ssm_parameter_t **parameters);
ssm_it_parameters_t *ssm_it_parameters_icdiff_new(ssm_parameter_t **parameters);

/* observed_template.c */
ssm_observed_t **_ssm_observed_new(int *observed_length);

/* diff_template.c */
void ssm_compute_diff(ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* ode_sde_template.c */
int ssm_step_ode(double t, const double X[], double f[], void *params);
void ssm_step_sde_no_dem_sto(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_step_sde_no_white_noise(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_step_sde_full(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_step_sde_no_dem_sto_no_white_noise(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* psr_template.c */
void ssm_psr_new(ssm_calc_t *calc);
void ssm_psr_free(ssm_calc_t *calc);
void ssm_step_psr(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* jac_template */
void ssm_eval_jac(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* Ht_template.c */
void ssm_eval_Ht(ssm_X_t *p_X, ssm_row_t *row, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* Q_template.c */
void ssm_eval_Q_no_dem_sto(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_eval_Q_no_env_sto(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_eval_Q_full(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);
void ssm_eval_Q_no_dem_sto_no_env_sto(const double X[], double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc);

/* step_ekf_template.c */
int ssm_step_ekf(double t, const double X[], double f[], void *params);



#endif
