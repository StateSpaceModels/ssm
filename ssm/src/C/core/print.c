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

void ssm_print_log(char *data)
{
    json_t *root;
    root = json_pack("{s,s,s,s}", "id", "log", "data", data);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
}

void ssm_print_warning(char *data)
{
    json_t *root;
    root = json_pack("{s,s,s,s}", "id", "wrn", "data", data);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
}

void ssm_print_err(char *data)
{
    json_t *root;
    root = json_pack("{s,s,s,s}", "id", "err", "data", data);
    json_dumpf(root, stderr, JSON_COMPACT); fprintf(stderr,"\n");
    fflush(stderr);
    json_decref(root);
}

void ssm_json_dumpf(FILE *stream, const char *id, json_t *data)
{
    json_t *root = json_pack("{s,s,s,o}", "id", id, "data", data);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
}

void ssm_print_X(FILE *stream, ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_row_t *row, const int index)
{
    //TODO handle t0 (ie before first data point)

    int i;

    ssm_state_t *state;
    ssm_observed_t *observed;
    double *X = p_X->proj;
    double t = (double) row->time;

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(index)); //j or m
    json_object_set_new(jout, "date", json_string(row->date));

    for(i=0; i<nav->states_sv_inc->length; i++){
        state = nav->states_sv_inc->p[i];
        json_object_set_new(jout, state->name, json_real(X[state->offset]));
    }

    for(i=0; i<nav->states_remainders->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(state->f_remainder(p_X, calc, t)));
    }

    for(i=0; i<nav->states_diff->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(state->f_inv(X[state->offset])));
    }

    for(i=0; i<nav->observed_length; i++){
        observed = nav->observed[i];
        json_object_set_new(jout, observed->name, json_real(observed->f_obs_mean(p_X, par, calc, t)));
    }

    char key[SSM_STR_BUFFSIZE];
    snprintf(key, SSM_STR_BUFFSIZE, "ran_%s", observed->name);
    for(i=0; i<nav->observed_length; i++){
        observed = nav->observed[i];
        json_object_set_new(jout, key, json_real(observed->f_obs_ran(p_X, par, calc, t)));
    }

    ssm_json_dumpf(stream, "X", jout);
}

/**
 * fitness is either log likelihood or sum of square
 */
void ssm_print_trace(FILE *stream, ssm_theta_t *theta, ssm_nav_t *nav, const double fitness, const int index)
{
    int i;
    ssm_parameter_t *parameter;

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(index)); // m

    for(i=0; i < nav->theta_all->length; i++) {
        parameter = nav->theta_all->p[i];
        json_object_set_new(jout, parameter->name, json_real(parameter->f_inv(gsl_vector_get(theta, i))));
    }

    json_object_set_new(jout, "fitness", isnan(fitness) ? json_null() : json_real(fitness));

    ssm_json_dumpf(stream, "trace", jout);
}



/**
 * computes standardized prediction residuals
 * res = (data-one_step_ahead_pred)/sqrt(var_one_step_ahead +var_obs)
 *
 * AND effective sample size.
 *
 * Note that this function is designed to be called only
 * ssm_data_t.length_nonan times by opposed to ssm_data_t.length (ie
 * when there is information)
 */
void ssm_print_pred_res(FILE *stream, ssm_X_t **J_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_row_t *row, ssm_fitness_t *fitness)
{
    int ts;
    double pred, var_obs, var_state, y, res;
    ssm_observed_t *observed;
    ssm_implementations_t implementation = nav->implementation;
    double t = (double) row->time;

    char key[SSM_STR_BUFFSIZE];

    //EKF specific
    ssm_X_t *X = *J_X;

    //non EKF
    int j;
    double kn, M2, delta, x;

    json_t *jout = json_object();
    json_object_set_new(jout, "date", json_string(row->date));

    for(ts=0; ts<row->ts_nonan_length; ts++) {

        observed = row->observed[ts];
        y = row->values[ts];

        if (implementation == SSM_EKF) {
            var_obs = observed->f_obs_var(X, par, calc, t);
            pred = observed->f_obs_mean(X, par, calc, t);
            var_state = observed->f_var_pred(X, par, calc, nav, t);
            res = (y - pred)/sqrt(var_state + var_obs);
        } else {
            kn=0.0;
            pred=0.0;
            var_obs=0.0;
            M2=0.0;

            for(j=0; j <fitness->J ; j++) {
                kn += 1.0;
                x = observed->f_obs_mean(J_X[j], par, calc, t);

                delta = x - pred;
                pred += delta/kn;
                M2 += delta*(x - pred);
                var_obs += observed->f_obs_var(J_X[j], par, calc, t);
            }

            var_state = M2/(kn - 1.0);
            var_obs /= ((double) fitness->J);

            res = (y - pred)/sqrt(var_state + var_obs);
        }

        snprintf(key, SSM_STR_BUFFSIZE, "pred_%s", observed->name);
        json_object_set_new(jout, key, json_real(pred));

        snprintf(key, SSM_STR_BUFFSIZE, "res_%s", observed->name);
        json_object_set_new(jout, key, (isnan(res)==1)? json_null(): json_real(res));
    }

    json_object_set_new(jout, "ess", json_real(fitness->ess_n));

    ssm_json_dumpf(stream, "predres", jout);
}


void ssm_print_hat(FILE *stream, ssm_hat_t *hat, ssm_nav_t *nav, ssm_row_t *row)
{
    int i;
    char key[SSM_STR_BUFFSIZE];

    ssm_state_t *state;
    ssm_observed_t *observed;

    json_t *jout = json_object();
    json_object_set_new(jout, "date", json_string(row->date));

    for(i=0; i< nav->states_sv_inc->length; i++) {
        state = nav->states_sv_inc->p[i];
        json_object_set_new(jout, state->name, json_real(hat->states[state->offset]));

        snprintf(key, SSM_STR_BUFFSIZE, "lower_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->states_95[state->offset][0]));

        snprintf(key, SSM_STR_BUFFSIZE, "upper_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->states_95[state->offset][1]));
    }

    for(i=0; i< nav->states_remainders->length; i++) {
        state = nav->states_remainders->p[i];
        json_object_set_new(jout, state->name, json_real(hat->remainders[state->offset]));

        snprintf(key, SSM_STR_BUFFSIZE, "lower_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->remainders_95[state->offset][0]));

        snprintf(key, SSM_STR_BUFFSIZE, "upper_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->remainders_95[state->offset][1]));
    }

    for(i=0; i< nav->states_diff->length; i++) {
        state = nav->states_diff->p[i];
        json_object_set_new(jout, state->name, json_real(hat->states[state->offset]));

        snprintf(key, SSM_STR_BUFFSIZE, "lower_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->states_95[state->offset][0]));

        snprintf(key, SSM_STR_BUFFSIZE, "upper_%s", state->name);
        json_object_set_new(jout, key, json_real(hat->states_95[state->offset][1]));
    }

    for(i=0; i< nav->observed_length; i++) {
        observed = nav->observed[i];

        json_object_set_new(jout, observed->name, json_real(hat->observed[observed->offset]));

        snprintf(key, SSM_STR_BUFFSIZE, "lower_%s", observed->name);
        json_object_set_new(jout, key, json_real(hat->observed_95[observed->offset][0]));

        snprintf(key, SSM_STR_BUFFSIZE, "upper_%s", observed->name);
        json_object_set_new(jout, key, json_real(hat->observed_95[observed->offset][1]));
    }

    ssm_json_dumpf(stream, "hat", jout);
}


/**
 * The key is to understand that: X_resampled[j] = X[select[j]] so select
 * give the index of the resample ancestor...
 * The ancestor of particle j is select[j]
 *
 * With n index: X[n+1][j] = X[n][select[n][j]]
 *
 * Other caveat: D_J_p_X are in [N_DATA+1] ([0] contains the initial conditions)
 * select is in [N_DATA], times is in [N_DATA]
 */
void ssm_sample_traj_print(FILE *stream, ssm_X_t ***D_J_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_data_t *data, ssm_fitness_t *fitness, const int index)
{
    int j_sel;
    int n, nn, indn;

    double ran, cum_weights;

    ssm_X_t *X_sel;

    ran=gsl_ran_flat(calc->randgsl, 0.0, 1.0);

    j_sel=0;
    cum_weights=fitness->weights[0];

    while (cum_weights < ran) {
        cum_weights += fitness->weights[++j_sel];
    }

    //print traj of ancestors of particle j_sel;

    //!!! we assume that the last data point contain information'
    X_sel = D_J_X[data->n_obs][j_sel]; // N_DATA-1 <=> data->indn_data_nonan[N_DATA_NONAN-1]
    ssm_print_X(stream, X_sel, par, nav, calc, data->rows[data->n_obs-1], index);

    //printing all ancesters up to previous observation time
    for(nn = (data->ind_nonan[data->n_obs_nonan-1]-1); nn > data->ind_nonan[data->n_obs_nonan-2]; nn--) {
        X_sel = D_J_X[ nn + 1 ][j_sel];
        ssm_print_X(stream, X_sel, par, nav, calc, data->rows[nn], index);
    }

    for(n = (data->n_obs_nonan-2); n >= 1; n--) {
        //indentifying index of the path that led to sampled particule
        indn = data->ind_nonan[n];
        j_sel = fitness->select[indn][j_sel];
        X_sel = D_J_X[ indn + 1 ][j_sel];

        ssm_print_X(stream, X_sel, par, nav, calc, data->rows[indn], index);

        //printing all ancesters up to previous observation time
        for(nn= (indn-1); nn > data->ind_nonan[n-1]; nn--) {
            X_sel = D_J_X[ nn + 1 ][j_sel];
            ssm_print_X(stream, X_sel, par, nav, calc, data->rows[nn], index);
        }
    }

    indn = data->ind_nonan[0];
    j_sel = fitness->select[indn][j_sel];
    X_sel = D_J_X[indn+1][j_sel];

    for(nn=indn; nn>=0; nn--) {
        X_sel = D_J_X[ nn + 1 ][j_sel];
        ssm_print_X(stream, X_sel, par, nav, calc, data->rows[nn], index);
    }

    //TODO nn=-1 (for initial conditions)

}


void ssm_print_ar(FILE *stream, ssm_adapt_t *adapt, const int index)
{

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(index)); // m

    json_object_set_new(jout, "ar", json_real(adapt->ar));
    json_object_set_new(jout, "ar_smoothed", json_real(adapt->ar_smoothed));
    json_object_set_new(jout, "eps", json_real(adapt->eps));

    ssm_json_dumpf(stream, "ar", jout);
}
