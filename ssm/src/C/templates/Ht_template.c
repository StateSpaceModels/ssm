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



void eval_Ht(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{

    double *X = p_X->proj;
    int i, j;
    gsl_matrix *Ft = calc->Ft;

    ssm_it_states_t *states_diff = nav->states_diff;
    ssm_it_states_t *states_inc = nav->states_inc;
    ssm_it_states_t *states_sv = nav->states_sv;


}
