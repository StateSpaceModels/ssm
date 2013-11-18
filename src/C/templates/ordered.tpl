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


{% for x, v in orders.items() %}
{% for o in v %}
#define ORDER_{{ o.name }} {{ o.order }}{% endfor %}{% endfor %}

#define E M_E
#define LN2 N_LN2
#define LN10 M_LN10
#define LOG2E M_LOG2E
#define LOG10E M_LOG10E
#define PI M_PI
#define SQRT1_2 M_SQRT1_2
#define SQRT2 M_SQRT2


{% block code %}
{% endblock %}

{% for x, v in orders.items() %}
{% for o in v %}
#undef ORDER_{{ o.name }}{% endfor %}{% endfor %}


#undef E
#undef LN2
#undef LN10
#undef LOG2E
#undef LOG10E
#undef PI
#undef SQRT1_2
#undef SQRT2

