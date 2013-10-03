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
 * Checking if initial conditions are valid
 */
ssm_err_code_t ssm_check_IC(ssm_par_t *par, ssm_calc_t *calc)
{

    cum_status = SSM_SUCCESS;

    /* checking remainders */
    {% for rem in parameters.remainders %}
    cum_status |= ({{ rem }} < 0) | SSM_ERROR_IC : cum_status;
    {% endfor %}

    /* checking state variables for populations without remainder */
    {% for ic in parameters.ics %}
    cum_status |= ({{ ic }} < 0) | SSM_ERROR_IC : cum_status;
    {% endfor %}

    return cum_status;
}
