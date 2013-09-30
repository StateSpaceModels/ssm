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

void ssm_compute_diff(ssm_X_t *p_X, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    {% if diff|length %}

    int i;
    ssm_it_states_t *it = nav->states_diff;

    double dt = p_X->dt;

    double _w[it->length];
    for(i=0; i<it->length; i++){
	_w[i] = gsl_ran_ugaussian(calc->randgsl);
    }

    {% for eq in diff %}
    p_X->proj[it->p[{{ loop.index0 }}]->offset] += sqrt(dt)*({{ eq }});{% endfor %}


    {% endif %}
}
