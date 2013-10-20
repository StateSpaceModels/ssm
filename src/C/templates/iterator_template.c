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

{% for t, its in iterators.items() %}
{% for k, orders in its.items() %}
ssm_it_{{ t }}s_t *ssm_it_{{ t }}s_{{ k }}_new(ssm_{{ t }}_t **{{ t }}s)
{
    ssm_it_{{ t }}s_t *it = _ssm_it_{{ t }}s_new({{ orders|length }});

    {% for x in orders %}
    it->p[{{ loop.index0 }}] = {{ t }}s[{{ x }}];{% endfor %}

    return it;
}
{% endfor %}{% endfor %}
