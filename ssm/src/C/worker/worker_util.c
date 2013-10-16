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

void ssm_zmq_send_par(void *socket, ssm_par_t *par, int zmq_options)
{   
    zmq_send(socket, par->data, par->size * sizeof (double), zmq_options);
}

void ssm_zmq_recv_par(ssm_par_t *par, void *socket)
{
    zmq_recv(socket, par->data, par->size * sizeof (double), 0);
}

void ssm_zmq_send_X(void *socket, ssm_X_t *X, int zmq_options)
{
    //dt
    zmq_send(socket, &(X->dt), sizeof (double), ZMQ_SNDMORE);
   
    //send proj
    zmq_send(socket, X->proj, X->length * sizeof (double), zmq_options);
}

void ssm_zmq_recv_X(ssm_X_t *X, void *socket)
{
    //dt
    zmq_recv(socket, &(X->dt), sizeof (double), 0);

    //proj
    zmq_recv(socket, X->proj, X->length * sizeof (double), 0);
}
