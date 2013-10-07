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



/**
 *  fill hat_95[2] with the 95% confidence interval (lower value in
 *  hat_95[0] and upper one in hat_95[1]). to_be_sorted is an array
 *  of the J particle to be sorted and weights their weights.
 *
 *  NOTE: if fitness is NULL, 1.0/calc->J is assumed as a weight
 */

void ssm_ci95(double *hat_95, ssm_calc_t *calc, ssm_fitness_t *fitness)
{
    int k;
    double weight_cum;

    int J = calc->J;
    const double *to_be_sorted = calc->to_be_sorted;
    size_t *index_sorted = calc->index_sorted;
    double *weights = (fitness) ? fitness->weights: NULL;

    //get the index of to_be_sorted.
    gsl_sort_index(index_sorted, to_be_sorted, 1, J); //to_be_sorted is not modified (i.e. not sorted in place), the index of the sorting are put in index_sorted.

    //cumulate sorted weight until we reach 2.5 % and take the corresponding value in to_be_sorted
    if(fitness){
	k=0; weight_cum = 0.0;
	while(weight_cum < 0.025) {
	    weight_cum += weights[index_sorted[k]];
	    k++;
	}
    } else {
	k = floor(0.025*J);
    }
    hat_95[0] = to_be_sorted[index_sorted[((k-1) <0) ? 0 : k-1]];

    //cumulate sorted weight until we reach 97.5 % and take the corresponding value in to_be_sorted
    if(fitness){
	k=0; weight_cum = 0.0;
	while(weight_cum < 0.975) {
	    weight_cum += weights[index_sorted[k]];
	    k++;
	} 
    } else {
	k = floor(0.975*J);
    }
    hat_95[1] = to_be_sorted[index_sorted[((k-1) <0) ? 0 : k-1]];
}



/**
 *  compute estimation of the 95% confidence interval (CI) the state variables
 *
 *  Note that the estimations are computed by a weighted average (each
 *  value is weighted by it's likelihood value). if fitness is NULL
 *  the weights are taken to be 1/calc->J
 *  is_J_par true J_par (ssm_par_t[J]) instead of &par
 */
void ssm_hat_eval(ssm_hat_t *hat, ssm_X_t **J_X, ssm_par_t **J_par, ssm_nav_t *nav, ssm_calc_t *calc, ssm_fitness_t *fitness, const double t, int is_J_par)
{

    int i, j;
    ssm_state_t *state;
    ssm_observed_t *observed;
    int offset;
    ssm_implementations_t implementation = nav->implementation;

    if (implementation == SSM_EKF) {
	int m = nav->states_sv->length + nav->states_inc->length + nav->states_diff->length;
	ssm_X_t *X = *J_X;
	ssm_par_t *par = *J_par;
	gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&X->proj[m], m, m);
	double rem, obs, var, grad;

	//sv and incidences
	for(i=0; i< nav->states_sv_inc->length; i++) {
	    offset = nav->states_sv_inc->p[i]->offset;
	    hat->states_95[offset][0] = X->proj[offset] - 1.96*gsl_matrix_get(&Ct.matrix,offset,offset);
	    hat->states_95[offset][1] = X->proj[offset] + 1.96*gsl_matrix_get(&Ct.matrix,offset,offset);
	}

	//remainders
	for(i=0; i< nav->states_remainders->length; i++) {
	    state = nav->states_remainders->p[i];
	    offset = state->offset;
	    rem = state->f_remainder(X, calc, t);
	    var = state->f_remainder_var(X, calc, nav, t);
	    hat->states_95[offset][0] = rem - 1.96*var;
	    hat->states_95[offset][1] = rem + 1.96*var;
	}

	//diffusions (we rely on a first-order Taylor approximation here)
	for(i=0; i<nav->states_diff->length; i++){
	    state = nav->states_diff->p[i];
	    offset = state->offset;
	    grad = state->f_inv_derivative(X->proj[offset]);
	    hat->states_95[offset][0] = state->f_inv(X->proj[offset]) - 1.96*pow(grad,2)*gsl_matrix_get(&Ct.matrix,offset,offset);
	    hat->states_95[offset][1] = state->f_inv(X->proj[offset]) + 1.96*pow(grad,2)*gsl_matrix_get(&Ct.matrix,offset,offset);
	}

	//observed
	for(i=0; i< nav->observed_length; i++) {
	    observed = nav->observed[i];
	    offset = observed->offset;
	    obs = observed->f_obs_mean(X, par, calc, t);
	    var = observed->f_obs_var(X, par, calc, t);
	    hat->states_95[offset][0] = obs - 1.96*var;
	    hat->states_95[offset][1] = obs + 1.96*var;
	}

    } else {
	//if fitness is NULL all the weights are set to 1.0/J
	int _zero = 0;
	double invJ[1] = {1.0 / (double) calc->J};
	int *j_par = (is_J_par) ? &j: &_zero;
	int *j_weights = (fitness) ? &j: &_zero;
	double *weights = (fitness) ? fitness->weights: invJ;
    
	//sv and incidences
	for(i=0; i< nav->states_sv_inc->length; i++) {
	    offset = nav->states_sv_inc->p[i]->offset;
	    hat->states[offset] = 0.0;
	    for(j=0; j<calc->J; j++) {
		calc->to_be_sorted[j] = J_X[j]->proj[offset]; //gsl_sort_index requires an array to be sorted and our particles are in J_X->proj so we use an helper array (calc->to_be_sorted)
		hat->states[offset] += calc->to_be_sorted[j]*weights[ *j_weights ];
	    }
	    ssm_ci95(hat->states_95[offset], calc, fitness);
	}
	
	//remainders
	for(i=0; i< nav->states_remainders->length; i++) {
	    state = nav->states_remainders->p[i];
	    offset = state->offset;
	    hat->remainders[offset] = 0.0;
	    for(j=0; j<calc->J; j++) {
		calc->to_be_sorted[j] = state->f_remainder(J_X[j], calc, t);
		hat->remainders[offset] += calc->to_be_sorted[j]*weights[ *j_weights ];
	    }
	    ssm_ci95(hat->remainders_95[offset], calc, fitness);
	}
	
	//diffusions
	for(i=0; i<nav->states_diff->length; i++){
	    state = nav->states_diff->p[i];
	    offset = state->offset;
	    hat->states[offset] = 0.0;
	    for(j=0; j<calc->J; j++) {
		calc->to_be_sorted[j] = state->f_inv(J_X[j]->proj[state->offset]);
		hat->states[offset] += calc->to_be_sorted[j]*weights[ *j_weights ];
	    }
	    ssm_ci95(hat->states_95[offset], calc, fitness);
	}
	
	//observed
	for(i=0; i< nav->observed_length; i++) {
	    observed = nav->observed[i];
	    offset = observed->offset;
	    hat->observed[offset] = 0.0;
	    for(j=0; j<calc->J; j++) {
		calc->to_be_sorted[j] = observed->f_obs_mean(J_X[j], J_par[ *j_par ], calc, t);
		hat->observed[offset] += calc->to_be_sorted[j]*weights[ *j_weights ];
	    }
	    ssm_ci95(hat->observed_95[offset], calc, fitness);
	}
	
    }
}
