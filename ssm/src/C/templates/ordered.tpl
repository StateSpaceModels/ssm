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

{% for o in orders.var %}
#define ORDER_{{ o }} {{ loop.index0 }}{% endfor %}

{% for o in orders.inc %}
#define ORDER_{{ o.name }} {{ o.order }}{% endfor %}

{% for o in orders.univ_rem %}
#define ORDER_{{ o.name }} {{ o.order }}{% endfor %}

{% for o in orders.diff %}
#define ORDER_{{ o }} {{ loop.index0 }}{% endfor %}

{% for o in orders.covariates %}
#define ORDER_{{ o }} {{ loop.index0 }}{% endfor %}


{% block code %}
{% endblock %}

{% for o in orders.var %}
#undef ORDER_{{ o }}{% endfor %}

{% for o in orders.inc %}
#undef ORDER_{{ o.name }}{% endfor %}

{% for o in orders.universe %}
#undef ORDER_{{ o.name }}{% endfor %}

{% for o in orders.diff %}
#undef ORDER_{{ o }}{% endfor %}

{% for o in orders.covariates %}
#undef ORDER_{{ o }}{% endfor %}

