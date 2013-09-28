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

unsigned int *ssm_u1_new(int n)
{
    int i;
    unsigned int *tab = malloc(n*sizeof (unsigned int));

    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i]=0;

    return tab;  
}

unsigned int **ssm_u2_new(int n, int p)
{
    int i;
    unsigned int **tab = malloc(n* sizeof (unsigned int *));

    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i] = ssm_u1_new(p);
  
    return tab;
}

void ssm_u2_free(unsigned int **tab, int n)
{
    int i;
    for(i=0; i<n; i++)
	free(tab[i]);

    free(tab);  
}

unsigned int ***ssm_u3_new(int n, int p1, int p2)
{
    int i;
    unsigned int ***tab = malloc(n* sizeof (unsigned int **));

    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i] = ssm_u2_new(p1, p2);

    return tab;  
}

void ssm_u3_free(unsigned int ***tab, int n, int p1)
{
    int i, j;

    for(i=0; i<n; i++)
	for(j=0; j<p1; j++)
	    free(tab[i][j]);

    for(i=0; i<n; i++)
	free(tab[i]);
    free(tab);
}

unsigned int ****ssm_u4_new(int n, int p1, int p2, int p3)
{
    int i;
    unsigned int ****tab = malloc(n* sizeof (unsigned int ***));

    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i] = ssm_u3_new(p1, p2, p3);

    return tab;
}

void ssm_u4_free(unsigned int ****tab, int n, int p1, int p2)
{
    int i, j, k;

    for(i=0; i<n; i++)
	for(j=0; j<p1; j++)
	    for(k=0; k<p2; k++)
		free(tab[i][j][k]);

    for(i=0; i<n; i++)
	for(j=0; j<p1; j++)
	    free(tab[i][j]);

    for(i=0; i<n; i++)
	free(tab[i]);

    free(tab);
}

unsigned int **ssm_u2_var_new(int n, unsigned int *p)
{
    int i;

    unsigned int **tab=malloc(n* sizeof (unsigned int *));
    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i] = ssm_u1_new(p[i]);

    return tab;
}


unsigned int ***ssm_u3_var_new(int n, unsigned int *p1, unsigned int **p2)
{
    int i;

    unsigned int ***tab = malloc(n* sizeof (unsigned int **));
    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
	tab[i]= ssm_u2u_var_set0(p1[i], p2[i]);

    return tab;
}


void ssm_u3_free(unsigned int ***tab, int n, unsigned int *p1)
{
    int i, j;

    for(i=0; i<n; i++)
	for(j=0; j<p1[i]; j++)
	    free(tab[i][j]);

    for(i=0; i<n; i++)
	free(tab[i]);

    free(tab);
}


unsigned int ***ssm_u3_varp1_new(int n, unsigned int *p1, int p2)
{
    int i,j;

    unsigned int ***tab = malloc(n* sizeof (unsigned int **));
    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
    {
	tab[i]= malloc(p1[i]* sizeof (unsigned int *));
	if(tab[i]==NULL)
	{
	    char str[STR_BUFFSIZE];
	    sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	    print_err(str);
	    exit(EXIT_FAILURE);
	}
    }

    for(i=0;i<n;i++)
	for(j=0;j<p1[i];j++)
	    tab[i][j]= ssm_u1_new(p2);

    return tab;
}


unsigned int ***ssm_u3_varp2_new(int n, unsigned int p1, unsigned int *p2)
{
    int i,j;

    unsigned int ***tab = malloc(n* sizeof (unsigned int **));
    if(tab==NULL)
    {
	char str[STR_BUFFSIZE];
	sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	print_err(str);
	exit(EXIT_FAILURE);
    }

    for(i=0;i<n;i++)
    {
	tab[i] = malloc(p1* sizeof (unsigned int *));
	if(tab[i]==NULL)
	{
	    char str[STR_BUFFSIZE];
	    sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);	
	    print_err(str);
	    exit(EXIT_FAILURE);
	}
    }

    for(i=0;i<n;i++)
	for(j=0;j<p1;j++)
	    tab[i][j]= ssm_u1_new(p2[i]);

    return tab;  
}
