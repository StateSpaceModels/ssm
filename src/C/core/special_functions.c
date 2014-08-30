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
 * The Heaviside step function or unit step function
 * x is typically t-t_intervention
 */
double heaviside(double x)
{    
    return (x < 0.0) ? 0.0 : 1.0;
}

/**
 * The ramp function 
 */
double ramp(double x)
{
    return (x >= 0) ? x : 0.0;
}

/**
 * Slowstep function : piecewise linear function equivalent to
 * integral of a rectangular function. Of the form
 * slowstep(x, d) = {   0     if x < 0
 *                      d     if x >= d
 *                      x     otherwise   }
 */
double slowstep(double x, double d)
{
    return ( x >= 0.0 ? ( x >= d ? d : x ) : 0.0 );
}

/**
 * Sigmoid function: decreases from 1 to 0.
 * x is typically t
 * shape controls the steep of the sigmoid (the greater the steeper)
 * shift controls how far from 0 the midpoint is shifted
*/
double sigmoid(double x, double shape, double shift) 
{
	return ( 1 / ( 1 + exp( - shape * ( x - shift ) ) ) );
}
