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
 * check if name is in it
 */
int ssm_in_par(ssm_it_parameters_t *it, const char *name)
{
    int i;
    for(i=0; i<it->length; i++){
	if(strcmp(it->p[i]->name, name) == 0){
	    return 1;
	}
    }

    return 0;
}
