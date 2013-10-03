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

char *ssm_c1_new(int n)
{
    char *tab = malloc(n*sizeof (char));

    if(tab==NULL)
    {
        char str[SSM_STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        ssm_print_err(str);
        exit(EXIT_FAILURE);
    }

    return tab;
}

char **ssm_c2_new(int n, int p)
{
    int i;
    char **tab = malloc(n* sizeof (char *));

    if(tab==NULL)
    {
        char str[SSM_STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        ssm_print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0 ; i<n ; i++)
        tab[i] = ssm_c1_new(p);

    return tab;
}

void ssm_c2_free(char **tab, int n)
{
    int i;
    for(i=0; i<n; i++)
        free(tab[i]);

    free(tab);
}
