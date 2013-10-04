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
    json_dumpf(root, stdout, 0); printf("\n");
    fflush(stdout);
    json_decref(root);
}

void ssm_print_warning(char *data)
{
    json_t *root;
    root = json_pack("{s,s,s,s}", "id", "wrn", "data", data);
    json_dumpf(root, stdout, 0); printf("\n");
    fflush(stdout);
    json_decref(root);
}

void ssm_print_err(char *data)
{
    json_t *root;
    root = json_pack("{s,s,s,s}", "id", "err", "data", data);
    json_dumpf(root, stderr, 0); fprintf(stderr,"\n");
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

void ssm_print_X(FILE *stream, ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_row_t *row, const int index, const double t)
{
    int i;

    ssm_state_t *state;
    ssm_observed_t *observed;
    double *X = p_X->proj;

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(index)); //j or m
    json_object_set_new(jout, "date", json_string(row->date));

    for(i=0; i<nav->states_sv->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(X[state->offset]));
    }

    for(i=0; i<nav->states_remainders->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(state->f_remainder(p_X, calc, t)));
    }

    for(i=0; i<nav->states_inc->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(X[state->offset]));
    }

    for(i=0; i<nav->states_diff->length; i++){
        state = nav->states_sv->p[i];
        json_object_set_new(jout, state->name, json_real(state->f_inv(X[state->offset])));
    }

    for(i=0; i<nav->observed_length; i++){
        observed = nav->observed[i];
        json_object_set_new(jout, observed->name, json_real(observed->f_obs_mean(p_X, par, calc, t)));
    }

    char prefix[SSM_STR_BUFFSIZE] = "ran_";
    for(i=0; i<nav->observed_length; i++){
        observed = nav->observed[i];
        json_object_set_new(jout, strcat(prefix, observed->name), json_real(observed->f_obs_ran(p_X, par, calc, t)));
    }

    ssm_json_dumpf(stream, "X", jout);
}


void ssm_print_trace(FILE *stream, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc, const int index, const double fitness)
{
    int i;
    ssm_parameter_t *parameter;

    json_t *jout = json_object();
    json_object_set_new(jout, "index", json_integer(index)); // m

    for(i=0; i < nav->theta_all->length; i++) {
        parameter = nav->theta_all->p[i];
        json_object_set_new(jout, parameter->name, json_real(parameter->f_par2user(gsl_vector_get(par, parameter->offset), par, calc)));
    }

    json_object_set_new(jout, "fitness", isnan(fitness) ? json_null() : json_real(fitness));

    ssm_json_dumpf(stream, "trace", jout);
}
