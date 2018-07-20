/*
METAPOPULATION MODEL
VIRAL REPLICATION HAPPENS IN MULTIPLE HOSTS WHO ARE CAPABLE OF TRANSMITTING THE VIRUS TO EACH OTHER IN A RANCOM NETWORK.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

//ran1
//poidev

int main()
{
	int i,j,k;
	int host_num = 1; // number of hosts in a system
	int kmax = 20; // max number of mutation in a segment. pop with kmax in a segment has fitness of 0.
	//pop: entire population including all metapops in each host. pop[0][1][2] = number of individual with 1 mutation in 1st segment and 2 mutation in 2nd segment in host 0.
	
	int*** pop = (int***) malloc(sizeof(int**)*host_num);
	for (i=0; i<kmax; i++)
	{
		pop[i] = (int**) malloc(sizeof(int*)*kmax);
		for (j=0; j<kmax; j++)
		{
			pop[i][j] = (int*) malloc(sizeof(int)*kmax);
		}
	}

	for (i=0; i<host_num; i++)
	{
		for (j=0; j<kmax; j++)
		{
			for (k=0; k<kmax; k++)
			{
				pop[i][j][k] = i+j+k;
				//printf("pop[%d][%d][%d]=%d\n",i,j,k,pop[i][j][k]);
			}
		}
	}





	return 0;
}