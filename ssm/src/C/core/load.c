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

void ssm_input2par(ssm_par_t *par, ssm_input_t *input, ssm_calc_t *calc, ssm_nav_t *nav)
{   
    ssm_it_parameters_t *it = nav->theta_all;
    
    int i;
    
    for(i=0; i< it->length; i++){
	gsl_vector_set(par, i, it->p[i]->f_user2par(gsl_vector_get(input, i), input, calc));
    }
}


void ssm_par2X(ssm_X_t *X, ssm_par_t *par, ssm_calc_t *calc, ssm_nav_t *nav)
{
    int i;

    ssm_it_states_t *sv = nav->states_sv;
    ssm_it_states_t *inc = nav->states_inc;
    ssm_it_states_t *diff = nav->states_diff;
    
    for(i=0; i<sv->length; i++){
	gsl_vector_set(X, sv->p[i]->offset, gsl_vector_get(par, sv->p[i]->ic->offset));
    }

    for(i=0; i<inc->length; i++){
	gsl_vector_set(X, inc->p[i]->offset, 0.0);
    }

    for(i=0; i<diff->length; i++){
	gsl_vector_set(X, diff->p[i]->offset, gsl_vector_get(par, diff->p[i]->ic->offset));
    }

}
