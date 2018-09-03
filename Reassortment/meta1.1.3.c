/*
METAPOPULATION MODEL
VIRAL REPLICATION HAPPENS IN MULTIPLE HOSTS WHO ARE CAPABLE OF TRANSMITTING THE VIRUS TO EACH OTHER IN A RANCOM NETWORK.
in 1.1.3, it uses simple one mutation accumulation equation.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

// function pre-decleration
//ran1
//poidev
void mutate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, double u, int kmax, int host_num, double* N2, double *N1, double factor[], int mutcap);
void reast(double**** pop2, int* curpop, int kmax, int host_num, double r, double* N2);
void repr(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double s, double* N2, double* N1, double* N, double c, double* K, long* seed);
float ran1(long *seed);
float gammln(float xx);
float poidev(float xm,long *idum);
float gasdev(long* idum);
double Nsum(int size,double N[]);
void record(double* N1, double* N2, double**** pop1, double**** pop2, int kmax, int host_num, int timestep, int krecord, int curpop1, int curpop2, int rep, int gen, FILE **fPointer);
void migrate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double tr, double mig, double* K, int K0, double* kf,int gen, double kf0, double kfsd);
void evol(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double evolq, double evolrate);
double fact(int num);
double poipmf(double l, int k);

#define PRINTF 0	
#define MUTPRINTF 0
#define REAPRINTF 0
#define REPPRINTF 0 
#define MIGPRINTF 0
#define EVOPRINTF 0
int main(int argc, char *argv[])
{
	// clock start
	//clock_t begin = clock();
	// prarmeter calling from os command
	char *destination = argv[1]; // directory to save data at
	char *timestep_s = argv[2]; // whether to record every generation
	char *krecord_s = argv[3]; // how to record k. (0=mean&min k, 1=min k, 2=histogram)
	char *untilext_s = argv[4]; // whether to run the simulation only until extinction (N = 0)
	char *rep_s = argv[5]; // number of repetition
	char *s_s = argv[6]; // selection coefficient
	char *N0_s = argv[7]; // base initial frequency
	char *K0_s = argv[8]; // carrying capacity in a host
	char *u_s = argv[9]; // mutation rate per segment
	char *gen_num_s = argv[10]; // number of generations
	char *c_s = argv[11]; // cost of multisegmentation
	char *r_s = argv[12]; // reassortment probability
	char *seed_s = argv[13]; // seed for random process
	char *host_num_s = argv[14]; // host number
	char *kmax_s = argv[15]; // maximum mutation amount in a segment
	char *pop2init_s = argv[16]; // information on the proportion of initial population to the N0 for 2 segments.
	char *pop2i_l_s = argv[17]; // subsidiary info for extracting pop2init info
	char *pop1init_s = argv[18]; // information on the proportion of initial population to the N0 for 1 segments.
	char *pop1i_l_s = argv[19]; // subsidiary info for extracting pop2init info
	char *tr_s = argv[20]; // transmission rate
	char *mig_s = argv[21]; // migration rate to the migration pool
	char *mutcap_s = argv[22]; // migration rate to the migration pool
	char *kf0_s = argv[23];
	char *kfsd_s = argv[24];
	char *evolrate_s = argv[25];
	char *evolq_s = argv[26];
	char *end1;

	int timestep = (int) strtol(timestep_s,&end1,10);
	int krecord = (int) strtol(krecord_s, &end1, 10);
	int untilext = (int) strtol(untilext_s,&end1,10);
	int rep = (int) strtol(rep_s,&end1,10);
	double s = (double) strtof(s_s, NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	int K0 = (int) strtol(K0_s,&end1,10);
	double u = (double) strtof(u_s,NULL);
	int gen_num = (int) strtol(gen_num_s,&end1,10);
	double c = (double) strtof(c_s,NULL);
	double r = (double) strtof(r_s,NULL);
	long seed = strtol(seed_s,&end1,10);
	int host_num = (int) strtol(host_num_s,&end1,10);
	int kmax = (int) strtol(kmax_s,&end1,10);
	int pop2i_l = (int) strtol(pop2i_l_s,&end1,10);
	int pop1i_l = (int) strtol(pop1i_l_s,&end1,10);
	double tr = (double) strtof(tr_s, NULL);
	double mig = (double) strtof(mig_s, NULL);
	int mutcap = (int) strtol(mutcap_s,&end1,10);
	double kf0 = (double) strtof(kf0_s, NULL);
	double kfsd = (double) strtof(kfsd_s, NULL);
	double evolrate = (double) strtof(evolrate_s, NULL);
	double evolq = (double) strtof(evolq_s, NULL);

	int k,i,j,m;
	//double pop1init[host_num];
	double pop2init[host_num];
	double pop1init[host_num];
	int len = 0;

	char* dubstr = (char*) malloc(sizeof(char)*10);
	for (i=0; i < pop2i_l; i++)
	{

		if (pop2init_s[i] != '~')
		{
			sprintf(dubstr,"%s%c",dubstr,pop2init_s[i]);
		}
		else
		{
			pop2init[len] = (double) strtof(dubstr, NULL);
			len++;
			sprintf(dubstr,"");
		}
	}
	len = 0;
	for (i=0; i < pop1i_l; i++)
	{
		if (pop1init_s[i] != '~')
		{
			sprintf(dubstr,"%s%c",dubstr,pop1init_s[i]);
		}
		else
		{
			pop1init[len] = (double) strtof(dubstr, NULL);
			len++;
			sprintf(dubstr,"");
		}
	}
	free(dubstr);
	printf("destination=%s, timestep=%d, krecord=%d, hostnum=%d, untilext=%d, kamx=%d, rep=%d, s=%.2f, N0=%d, K0=%d, u=%.5f, gen_num=%d, c=%.2f, r=%.2f, tr=%.5f, mig=%.5f, mutcap=%d, seed=%lu, kf=%.2f, kfsd=%.2f, evolrate=%.2f, evolq=%.2f\n",destination,timestep,krecord,host_num,untilext,kmax,rep,s,N0,K0,u,gen_num,c,r,tr,mig,mutcap,seed,kf0,kfsd,evolrate,evolq);

	//check if the destination folder exists and if not, make one.
	char* dest2 = (char*) malloc(sizeof(char)*50);
	sprintf(dest2,"./data/%s",destination);
	struct stat st = {0};
	if (stat(dest2,&st) == -1)
	{
		mkdir(dest2, 0700);
	}

	char* filename = (char*) malloc(sizeof(char)*1000);
	sprintf(filename,"%s/m1.1.3s_%d,%d,%d,%.3f,%d,%d,%.5f,%d,%.2f,%.2f,%d,%d,%.5f,%.5f,%.2f,%.2f,%.2f,%.2f(0).csv",dest2,timestep,krecord,rep,s,N0,K0,u,gen_num,c,r,kmax,host_num,mig,tr,kf0,kfsd,evolrate,evolq);
	int filenum = 0;

	while ( access(filename, F_OK) != -1) 
	{
		filenum += 1;
		sprintf(filename,"%s/m1.1.3s_%d,%d,%d,%.3f,%d,%d,%.5f,%d,%.2f,%.2f,%d,%d,%.5f,%.5f,%.2f,%.2f,%.2f,%.2f(%d).csv",dest2,timestep,krecord,rep,s,N0,K0,u,gen_num,c,r,kmax,host_num,mig,tr,kf0,kfsd,evolrate,evolq,filenum);
	}
	FILE* fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);
	char* str;
	if (krecord == 1 || krecord == 0)
	{
		str = (char*) malloc(sizeof(char)*(20*(host_num+1)*2) + 100); // (words + comma)*(host_num)*(seg1 and 2) + safety buffer
	}
	else if (krecord == 2)
	{
		str = (char*) malloc(sizeof(char)*(20*(2*kmax+1)*2) + 100); // (words + comma)*(host_num)*(seg1 and 2) + safety buffer
	}
	sprintf(str,"");
	if(krecord != 2)
	{
		for(i=0; i<=host_num; i++)
		{
			sprintf(str,"%s,pop1.%d,pop2.%d,k1.%d,k2.%d",str,i,i,i,i);	
		}
	}
	else
	{
		sprintf(str,"%s,pop1.0,pop2.0",str);
		for (i=0; i<=kmax*2; i++)
		{
			sprintf(str,"%s,p1k%d,p2k%d",str,i,i);
		}
	}
	if (timestep)
	{
		fprintf(fPointer,"rep,gen%s\n",str);
		//printf("rep,gen%s\n",str);
	}
	else
	{
		fprintf(fPointer,"rep%s\n",str);
		//printf("rep,gen%s\n",str);
	}
	free(str);
	// todo:
	
	// set the parameters with os command.
	// set timer, file
	
	int gen,repe; //current generation and repetition.
	int curpop1, curpop2; //current population index ur working witih.
	double**** pop2; //population of 2 segments
	double**** pop1; //population of 1 segments
	double* N = (double*) malloc(sizeof(double)*(host_num+1));
	double* N2 = (double*) malloc(sizeof(double)*(host_num+1)); // current populationop size of 2segs. N2[i] is a pop size of host i. N2[0] is the total population size of 2segs.
	double* N1 = (double*) malloc(sizeof(double)*(host_num+1)); // current populationop size of 1segs. N1[i] is a pop size of host i. N1[0] is the total population size of 1segs.
	double* K = (double*) malloc(sizeof(double)*(host_num));
	double* kf = (double*) malloc(sizeof(double)*(host_num));
	for (i=0; i<host_num; i++)
	{
		kf[i] = gasdev(&seed)*kfsd + kf0;
	}
	double factor[2*kmax]; // probability of getting n mutations organized in an array
	//double record; // record of mutation amount in a pop of a host (mean or minimum depending on krecord)
	//pop2: entire population of 2 segments including all metapops in each host. pop[i][1][1][2] = number of individual with 1 mutation in 1st segment and 2 mutation in 2nd segment in host 0.
	//pop1: same thing as pop2 but of 1 segments. pop[i][1][1] = number of individual 1 with 1 mutation in its genome.
	// the first index of pop indicates where the current population is stored. first index alternates between 0 and 1 in order to process the generation step.
	//host 0 (pop1[i][0][j] or pop2[i][0][j][k]) is a migration pool.
	pop2 = (double****) malloc(sizeof(double***)*2);
	for (m=0; m<2; m++)
	{
		pop2[m] = (double***) malloc(sizeof(double**)*(host_num+1));
		for (i=0; i<=host_num; i++)
		{
			pop2[m][i] = (double**) malloc(sizeof(double*)*(kmax+1));
			for (j=0; j<=kmax; j++)
			{
				pop2[m][i][j] = (double*) malloc(sizeof(double)*(kmax+1));
			}
		}
	}


	pop1 = (double****) malloc(sizeof(double***)*2);
	for (m=0; m<2; m++)
	{
		pop1[m] = (double***) malloc(sizeof(double**)*(host_num+1));
		for (i=0; i<=host_num; i++)
		{
			pop1[m][i] = (double**) malloc(sizeof(double*)*(kmax+1));
			for (j=0; j<=kmax; j++)
			{
				pop1[m][i][j] = (double*) malloc(sizeof(double)*(kmax+1));
			}
		}
	}

	for (i=0; i<=2*kmax; i++){
		factor[i] = poipmf(2*u,i); // probability of choosing i amount of mutation in a generation step.
	}

	double mutate_time = 0;
	double reast_time = 0;
	double repr_time = 0;
	double migr_time = 0;
	double reco_time = 0;
	clock_t begin, end;

	for (repe=0; repe < rep; repe++)
	{
		if (rep > 100)
		{
			if (repe % 100 == 0)
			{
				printf("\rREP = %d",repe);
				fflush(stdout);
			}
		}
		else
		{
			printf("\rREP = %d",repe);
			fflush(stdout);
		}
		for(i=0; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					{
						pop2[0][i][j][k] = 0;
						pop1[0][i][j][k] = 0;
					}
				}
			}
			if (i!=0)
			{
				pop1[0][i][0][0] = (double) N0*pop1init[i-1];
				N1[i] = (double) N0*pop1init[i-1];
				if (evolrate > 0)
				{
					N2[i] = 0;
				}
				else
				{
					pop2[0][i][0][0] = (double) N0*pop2init[i-1]; //N0 of virus with 0 mutations at initial condition.
					N2[i] = (double) N0*pop2init[i-1];
				}
				N[i] += N1[i] + N2[i];
				K[i-1] = (double) K0;
			}
		}

		N1[0] = Nsum(host_num,N1);
		N2[0] = Nsum(host_num,N2);
		N[0] = N1[0] + N2[0];
		for(i=0; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					{
						//printf("pop2[0][%d][%d][%d]=%.2f\n",i,j,k,pop2[0][i][j][k]);
					}
				}
			}
		}

		//N1[0] = Nsum(host_num,N1);
		curpop1 = 0;
		curpop2 = 0;
		for (gen=0; gen < gen_num; gen++)
		{
			if ( N2[0] > 0 || N1[0] > 0)
			{
				begin = clock();
				mutate(pop2, pop1, &curpop2, &curpop1, u, kmax, host_num, N2, N1, factor, mutcap);
				end = clock();
				mutate_time += (double) (end - begin)/ CLOCKS_PER_SEC;
				begin = clock();
				reast(pop2, &curpop2, kmax, host_num, r, N2);
				end = clock();
				reast_time += (double) (end - begin)/ CLOCKS_PER_SEC;
				begin = clock();
				repr(pop2, pop1, &curpop2, &curpop1, kmax, host_num,s,N2,N1,N,c,K,&seed);
				end = clock();
				repr_time += (double) (end - begin)/ CLOCKS_PER_SEC;				
				begin = clock();
				migrate(pop2, pop1, &curpop2, &curpop1, kmax, host_num, N2, N1, N, &seed, tr, mig, K, K0, kf, gen, kf0, kfsd);
				end = clock();
				migr_time += (double) (end - begin)/ CLOCKS_PER_SEC;
				evol(pop2,pop1, &curpop2, &curpop1, kmax, host_num, N2, N1, N, &seed, evolq, evolrate);
				N[0] = Nsum(host_num,N);
				//printf("N[0]=%.3f\n",N[0]);
				//printf("N=%.3f N[1]=%.3f, N[2]=%.3f N[3]=%.3f N[4]=%.3f N[5]=%.3f\n",N[0],N[1],N[2],N[3],N[4],N[5]);
			}
			else if (N2[0] == 0 || N1[0] == 0)
			{
				if (untilext == 1)
				{
					//printf("break activated at gen=%d\n",gen);
					break;
				}
			}

			if(PRINTF)
			{
				float count1, count2;
				count1 = 0;
				count2 = 0;
				printf("after repr outside the function\n");
				for (i=1; i<=host_num; i++)
				{
					for (j=0; j<=kmax; j++)
					{
						for (k=0; k<=kmax; k++)
						{
							if (pop2[curpop2][i][j][k]>0)
							{
								printf("pop2[%d][%d][%d][%d]=%.3f\n",curpop2,i,j,k,pop2[curpop2][i][j][k]);
							}
							count2 += pop2[curpop2][i][j][k];
						}
					}
					for (j=0; j<=kmax; j++)
					{
						for (k=0; k<=kmax; k++)
						{
							if (pop1[curpop1][i][j][k]>0)
							{
								printf("pop2[%d][%d][%d][%d]=%.3f\n",curpop2,i,j,k,pop1[curpop1][i][j][k]);
							}
							count1 += pop1[curpop1][i][j][k];
						}
					}
				}
				printf("count2=%.3f\n",count2);
				printf("count1=%.3f\n",count1);
			}
			
			//printf("gen=%d\n",gen);
			//printf("-----------------\n");
			// record it to the file
			if(timestep == 1)
			{
				begin = clock();
				record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
				end = clock();
				reco_time += (double) (end - begin)/ CLOCKS_PER_SEC;				
			}

		}
		if(timestep == 0)
		{
			begin = clock();
			record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
			end = clock();
			reco_time += (double) (end - begin)/ CLOCKS_PER_SEC;				
		}
		// erase pop data to start over simulation
		for (m=0; m<2; m++)
		{
			for (i=0; i<=host_num; i++)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						pop2[m][i][j][k] = 0;
						pop1[m][i][j][k] = 0;
					}
				}
			}
		}
	}
	// freeing pop2
	for (m=0; m<2; m++)
	{
		for (i=0; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				free(pop2[m][i][j]);
				free(pop1[m][i][j]);
			}
			free(pop2[m][i]);
			free(pop1[m][i]);
			}
		free(pop2[m]);	
		free(pop1[m]);	
	}
	free(pop2);
	free(pop1);
	free(N2);
	free(N1);
	free(N);
	// make mutation, recombination, reproduction into a single process.
	fclose(fPointer);
	// clock_t end = clock();
	//double time_spent = (end - begin)/ CLOCKS_PER_SEC;
	printf("\n");
	printf("mutate=%.3f, reast=%.3f, repr=%.3f, migr=%.3f\n",mutate_time, reast_time, repr_time, migr_time);
	//printf("time spend was %.2f minutes\n", time_spent/60.0);
	return 0;
}


void mutate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, double u, int kmax, int host_num, double* N2, double *N1, double factor[], int mutcap)
{
	int s2m,s2m2, s1m, s1m2, i,j,k,l,l2,l3;
	if (*curpop2 == 0)
	{
		s2m2 = *curpop2;
		s2m = 1;
		*curpop2 = s2m;
	}
	else 
	{
		s2m2 = *curpop2;
		s2m = 0;
		*curpop2 = s2m;
	}

	if (*curpop1 == 0)
	{
		s1m2 = *curpop1;
		s1m = 1;
		*curpop1 = s1m;
	}
	else 
	{
		s1m2 = *curpop1;
		s1m = 0;
		*curpop1 = s1m;
	}

	double count1, count2;

	if(MUTPRINTF)
	{
		count2 = 0;
		count1 = 0;
		printf("before mutation inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m2][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
					}
					count2 += pop2[s2m2][i][j][k];
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m2][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
					}
					count1 += pop1[s1m2][i][j][k];
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
		printf("\n");
	}	

	double factor2 = 0;
	int left; // max of number of mutations that n(l,k) can give rise to.
	int select;
	int cap;
	for (i=1; i<=host_num; i++)
	{
		if (N2[i] > 0) // 2 seg mutation process
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// going to sum all the mutation_rate*n(l,k)
					left = 2*kmax - (j + k);
					pop2[s2m][i][j][k] += pop2[s2m2][i][j][k];
					if (left < mutcap) // if possible number of accumulation is smaller than mutcap, the max number of accumulation is the cap.
					{
						cap = left; 
					}
					else
					{
						cap = mutcap;
					}
					for(l=1; l<=cap; l++)
					{
						factor2 = factor[l]*pop2[s2m2][i][j][k];
						pop2[s2m][i][j][k] -= factor2;
						for(l2=0; l2<=l; l2++)
						{
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								if(l <= kmax - k && l <= kmax - j)
								{
									pop2[s2m][i][j+l2][k+l3] += factor2/(l + 1);
								}
								else if(l <= kmax - k || l <= kmax - j)
								{
									if (k > j)
									{
										select = k;
									}
									else
									{
										select = j;
									}
									pop2[s2m][i][j+l2][k+l3] += factor2/(kmax - select + 1);
								}
								else
								{
									pop2[s2m][i][j+l2][k+l3] += factor2/(2*kmax - k - j - l + 1);
								}
							}
						}
					}
				}
			}
		}
		if (N1[i] > 0) // 2 seg mutation process
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// going to sum all the mutation_rate*n(l,k)
					left = 2*kmax - (j + k);
					pop1[s1m][i][j][k] += pop1[s1m2][i][j][k];
					if (left < mutcap) // if possible number of accumulation is smaller than mutcap, the max number of accumulation is the cap.
					{
						cap = left; 
					}
					else
					{
						cap = mutcap;
					}
					for(l=1; l<=cap; l++)
					{
						factor2 = factor[l]*pop1[s1m2][i][j][k];
						pop1[s1m][i][j][k] -= factor2;
						for(l2=0; l2<=l; l2++)
						{
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								if(l <= kmax - k && l <= kmax - j)
								{
									pop1[s1m][i][j+l2][k+l3] += factor2/(l + 1);
								}
								else if(l <= kmax - k || l <= kmax - j)
								{
									if (k > j)
									{
										select = k;
									}
									else
									{
										select = j;
									}
									pop1[s1m][i][j+l2][k+l3] += factor2/(kmax - select + 1);
								}
								else
								{
									pop1[s1m][i][j+l2][k+l3] += factor2/(2*kmax - k - j - l + 1);
								}
							}
						}
					}
				}
			}
		}
	}

	if(MUTPRINTF)
	{
		count2 = 0;
		count1 = 0;
		printf("after mutation inside the function\n");
		printf("\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{	
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
					}
					//pop2[m2][i][j][k] = 0;
					count2 += pop2[s2m][i][j][k];
				}
			}
			for (j=0; j<=kmax; j++)
			{	
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
					}
					//pop2[m2][i][j][k] = 0;
					count1 += pop1[s1m][i][j][k];
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}
}


void reast(double**** pop2, int* curpop, int kmax, int host_num, double r, double* N2)
{
	int m,m2,i,j,k;
	int jj;
	if (*curpop == 0)
	{
		m2 = *curpop;
		m = 1;
		*curpop = m;
	}
	else 
	{
		m2 = *curpop;
		m = 0;
		*curpop = m;
	}
	double count2;
	double mean2 = 0;
	if(REAPRINTF)
	{
		count2 = 0;
		printf("before reast inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					//printf("pop2[%d][%d][%d][%d]=%.2f\n",m2,i,j,k,pop2[m2][i][j][k]);
					count2 += pop2[m2][i][j][k];
					mean2 += pop2[m2][i][j][k]*(j+k);
				}
			}
		}
		mean2 = mean2 / count2;
		printf("count2=%.3f\n",count2);
		printf("mean2=%.3f\n",mean2);

	}	
	double jp[kmax+1];
	double kp[kmax+1];
	for (i=1; i<=host_num; i++)
	{
		if (N2[i] > 0)
		{
			// calculate all the proportions first
			for(j=0; j<=kmax; j++)
			{
				kp[j] = 0;
				jp[j] = 0;
				for(jj=0; jj<=kmax; jj++)
					{
						kp[j] += pop2[m2][i][jj][j];
						jp[j] += pop2[m2][i][j][jj];
					}
					kp[j] = kp[j]/N2[i];
					jp[j] = jp[j]/N2[i];
			}
			// those that didn't recombine remain, and the recombined ones get added.
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop2[m][i][j][k] = pop2[m2][i][j][k]*(1 - r);
					pop2[m][i][j][k] += N2[i]*kp[k]*jp[j]*r;
				}
			}
		}
	}


	if (REAPRINTF)
	{
		count2 = 0;
		mean2 = 0;
		printf("after reast inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop2[m2][i][j][k] = 0;
					//printf("pop2[%d][%d][%d][%d]=%.2f\n",m,i,j,k,pop2[m][i][j][k]);
					count2 += pop2[m][i][j][k];
					mean2 += pop2[m][i][j][k]*(j+k);
				}
			}
		}
		mean2 = mean2/count2;
		printf("count2=%.3f\n",count2);
		printf("mean2=%.3f\n",mean2);
	}
}

void repr(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double s, double* N2, double* N1, double* N, double c, double* K, long* seed)
{
	int s2m,s2m2,s1m,s1m2,i,j,k;
	if (*curpop2 == 0)
	{
		s2m2 = *curpop2;
		s2m = 1;
		*curpop2 = s2m;
	}
	else 
	{
		s2m2 = *curpop2;
		s2m = 0;
		*curpop2 = s2m;
	}
	if (*curpop1 == 0)
	{
		s1m2 = *curpop1;
		s1m = 1;
		*curpop1 = s1m;
	}
	else 
	{
		s1m2 = *curpop1;
		s1m = 0;
		*curpop1 = s1m;
	}

	//double newN2, newN1;
	double count1, count2;
	if (REPPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("before repr inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m2][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
						count2 += pop2[s2m2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m2][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
						count1 += pop1[s1m2][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}
	double poirate;
	double g = 1;
	for (i=1; i<=host_num; i++)
	{
		if(N2[i] > 0)
		{
			//newN2 = 0; //reset N2 to 0
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// if total mutation number = 2*kmax the reproduction rate is 0.
					if (k + j == 2*kmax)
					{
						pop2[s2m][i][j][k] = poidev(0,seed);

					}
					else
					{
						poirate = pop2[s2m2][i][j][k] * pow((1 - s), (k + j)) * (1 - c) * ((double)(1+g) / (1.0 + g*(N[i]/K[i-1])));
						pop2[s2m][i][j][k] = poidev(poirate,seed);
					}
				}
			}
		}
		if(N1[i] > 0)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// if total mutation number = 2*kmax the reproduction rate is 0.
					if (k + j == 2*kmax)
					{
						pop1[s1m][i][j][k] = poidev(0,seed);
					}
					else
					{
						poirate = pop1[s1m2][i][j][k] * pow((1 - s), (k + j)) * (1 - c) * ((double)2 / (1.0 + (N[i]/K[i-1])));
						pop1[s1m][i][j][k] = poidev(poirate,seed);
					}
				}
			}
		}
	}
	if(REPPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("after repr inside the function\n");
	}

	// erase m2 info in order to make mutate function work on the next generation.
	for (i=0; i<=host_num; i++)
	{
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				//pop2[s2m2][i][j][k] = 0;
				//if (pop2[m][i][j][k] >0 )
				{
					if (REPPRINTF)
					{
						if (pop2[s2m][i][j][k] > 0)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
							count2 += pop2[s2m][i][j][k];
						}
					}
				}
			}
		}
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				if (REPPRINTF)
				{
					if (pop1[s1m][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
						count1 += pop1[s1m][i][j][k];
					}
				}
			}
		}

	}
	if (REPPRINTF)
	{
		printf("count1=%.3f\n",count1);
		printf("count2=%.3f\n",count2);
	}
}

void migrate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double tr, double mig, double* K, int K0, double* kf, int gen, double kf0, double kfsd)
{
	int s2m,s2m2,s1m,s1m2,i,j,k;
	double count1, count2;
	if (*curpop2 == 0)
	{
		s2m2 = *curpop2;
		s2m = 1;
		*curpop2 = s2m;
	}
	else 
	{
		s2m2 = *curpop2;
		s2m = 0;
		*curpop2 = s2m;
	}
	if (*curpop1 == 0)
	{
		s1m2 = *curpop1;
		s1m = 1;
		*curpop1 = s1m;
	}
	else 
	{
		s1m2 = *curpop1;
		s1m = 0;
		*curpop1 = s1m;
	}
	
	double newN2, newN1;

	if (MIGPRINTF)
	{
		count1 = 0;
		count2 = 0;
		for (i=0; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m2][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
						count2 += pop2[s2m2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m2][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
						count1 += pop1[s1m2][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}

	// going into the migration pool (pop[m][0])
	for(i=1; i<= host_num; i++)
	{
		if (N2[i] > 0)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop2[s2m][i][j][k] = pop2[s2m2][i][j][k];
					pop2[s2m][i][j][k] -= pop2[s2m2][i][j][k] * mig;
					pop2[s2m][0][j][k] += pop2[s2m2][i][j][k] * mig;
					//printf("pop2[s2m2][%d][%d][%d]=%.2f,mig=%.2f\n",i,j,k,pop2[s2m2][i][j][k],mig);
				}
			}
		}
		if (N1[i] > 0)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop1[s1m][i][j][k] = pop1[s1m2][i][j][k];
					pop1[s1m][i][j][k] -= pop1[s1m2][i][j][k] * mig;
					pop1[s1m][0][j][k] += pop1[s1m2][i][j][k] * mig;
				}
			}
		}
	}

	// transmission to hosts from migration pool
	for(i=1; i<= host_num; i++)
	{
		newN2 = 0; //reset N2 to 0
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				pop2[s2m][i][j][k] += poidev(pop2[s2m][0][j][k]/host_num*tr,seed);
				newN2 += pop2[s2m][i][j][k];

			}
		}
		N2[i] = newN2;
		newN1 = 0;
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				pop1[s1m][i][j][k] += poidev(pop1[s1m][0][j][k]/host_num*tr,seed);
				newN1 += pop1[s1m][i][j][k];

			}
		}
		N1[i] = newN1;
		N[i] = N2[i] + N1[i];
		if (K[i-1] > 0.00000001)
		{
			K[i-1] -= kf[i];
			//printf("kf=%.3f\n",kf);
			//printf("K[%d]inside=%.2f\n",i-1,K[i-1]);
		}
		if (K[i-1] <= 0)
		{
			K[i-1] = 0.00000001;
		}
		if (N[i] == 0)
		{
			kf[i-1] = gasdev(seed)*kfsd + kf0;
			K[i-1] = (double) K0;
		}
	}
	//double Kmean=0, Kvar=0;
	//int Ncount = 0;
	//if (gen == 200 || gen == 300 || gen == 400)
	/*{
		for (i=0; i< host_num; i++)
		{
			/
			if (N[i] > 0)
			{
				Ncount++;
			}
			if (gen == 200 || gen == 300 || gen == 400)
			{
				printf("K[%d]=%.0f\n",i,K[i]);
			}
			//Kmean += K[i];
		}
		if (Ncount == host_num && gen < 70)
		{
			printf("all hosts infected at gen=%d\n",gen);
		}

		Kmean = Kmean/host_num;
		for (i=0; i< host_num; i++)
		{
			Kvar += pow((K[i] - Kmean),2);
		}
		Kvar = Kvar/host_num;
		printf("Kmean = %.2f, Kvar = %.2f\n",Kmean,Kvar);
	}*/


	N2[0] = Nsum(host_num,N2);
	N1[0] = Nsum(host_num,N1);
	
	// erase migration info

	if(MIGPRINTF)
	{
		printf("migration pool\n");
	}
	for (j=0; j<=kmax; j++)
	{
		for (k=0; k<=kmax; k++)
		{
			if(MIGPRINTF)
			{
				if(pop1[s1m][0][j][k]>0)
				{
					printf("pop1[%d][0][%d][%d]=%.2f\n",s1m,j,k,pop1[s1m][0][j][k]);
				}
			}
			pop1[s1m][0][j][k] = 0;
		}
	}
	for (j=0; j<=kmax; j++)
	{
		for (k=0; k<=kmax; k++)
		{
			if(MIGPRINTF)
			{
				if(pop2[s2m][0][j][k]>0)
				{
					printf("pop2[%d][0][%d][%d]=%.2f\n",s2m,j,k,pop2[s2m][0][j][k]);
				}
			}
			pop2[s2m][0][j][k] = 0;
		}
	}

	if(MIGPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("after migr inside the function\n");
	}
	// erase m2 info in order to make mutate function work on the next generation.
	for (i=1; i<=host_num; i++)
	{
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				pop2[s2m2][i][j][k] = 0;
				pop1[s1m2][i][j][k] = 0;
				//if (pop2[m][i][j][k] >0 )
				{
					if (MIGPRINTF)
					{
						if(pop2[s2m][i][j][k] > 0)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
							count2 += pop2[s2m][i][j][k];
						}
						if(pop1[s1m][i][j][k] > 0)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
							count1 += pop1[s1m][i][j][k];
						}
					}
				}
			}
		}
	}
	if(MIGPRINTF)
	{
		printf("count1=%.2f\n",count1);
		printf("count2=%.2f\n",count2);
	}
}

void evol(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double evolq, double evolrate)
{
	// the class with the biggest size in pop1 transfers to pop2 with the proportion of evolq
	int s2,s1,i,j,k,l;
	double count1, count2;
	s2 = *curpop2;
	s1 = *curpop1;
	
	if (EVOPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("before evol inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",s2,i,j,k,pop2[s2][i][j][k]);
						count2 += pop2[s2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1,i,j,k,pop1[s1][i][j][k]);
						count1 += pop1[s1][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}

	double max = 0;
	int ii, jj, kk;
	for (l=0; l<(int)poidev(evolrate,seed); l++)
	{	
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1][i][j][k] > max)
					{
						ii = i;
						jj = j;
						kk = k;
						max = pop1[s1][i][j][k];
					}
				}
			}
		}
		pop2[s2][ii][jj][kk] = pop1[s1][ii][jj][kk]*evolq;
		pop1[s1][ii][jj][kk] -= pop1[s1][ii][jj][kk]*evolq;
	}
	if (EVOPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("after evol inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2][i][j][k] > 0)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",s2,i,j,k,pop2[s2][i][j][k]);
						count2 += pop2[s2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1][i][j][k] > 0)
					{
						printf("pop1[%d][%d][%d][%d]=%.2f\n",s1,i,j,k,pop1[s1][i][j][k]);
						count1 += pop1[s1][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}

}

double Nsum(int size,double N[])
{
  int i;
  double sum = 0;
  for(i=1; i<=size; i++){
    sum += N[i];
  }
  return sum;
}

void record(double* N1, double* N2, double**** pop1, double**** pop2, int kmax, int host_num, int timestep, int krecord, int curpop1, int curpop2, int rep, int gen, FILE **fPointer)
{
	int i,j,k;
	double krecord1, krecord2;
	double krecord1t, krecord2t;
	double p1hist[2*kmax+1], p2hist[2*kmax+1];
	for (j=0; j<=2*kmax; j++) // pop1
	{
		p2hist[j] = 0;
		p1hist[j] = 0;
	}

	char* str;
	if (krecord == 0 || krecord == 1)
	{
		str = (char*) malloc(sizeof(char)*(10+10+2*2*20*host_num));
	}
	else if (krecord == 2)
	{
		str = (char*) malloc(sizeof(char)*((2*kmax+1)*2*20 + 50)); //(number of mutations)*(1&2seg)*(digits for pop size + commas) + buffer
 	}
	sprintf(str,"");
	//char* str = (char*) malloc(sizeof(char)*(10+10+2*2*20*host_num)); // repe_num + gen_num + (1&2seg)*(k&pop)*(value)*(host_num)
	if(krecord == 0) // record mean
	{
		krecord1t = 0; // record of entire population for 1 seg
		krecord2t = 0; // `` for 2 seg

		for (i=1; i<=host_num; i++)
		{
			krecord1 = 0; // record for 1 host for 1seg
			krecord2 = 0; // for 2 seg
			if (N2[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						krecord2 += pop2[curpop2][i][j][k]/N2[i] * (j + k);
					}
				}
				krecord2t += krecord2*N2[i]/N2[0];
			}
			else
			{
				krecord2 = -1.0;
			}
			if (N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						krecord1 += pop1[curpop1][i][j][k]/N1[i] * (j + k);
					}
				}
				krecord1t += krecord1*N1[i]/N1[0];
			}
			else
			{
				krecord1 = -1.0;
			}
			sprintf(str,"%s,%.2f,%.2f,%.2f,%.2f",str,N1[i],N2[i],krecord1,krecord2);
		}
	}
	else if(krecord == 1)// record minimum
	{
		krecord2 = kmax*2 + 1;
		krecord1 = kmax*2 + 1;
		for (i=1; i<=host_num; i++)
		{
			if (N2[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if((j + k) < krecord2)
						{
							if(pop2[curpop2][i][j][k] > 0)
							{
								krecord2 = j + k;
							}
						}
					}
				}
				if (krecord2t > krecord2)
				{
					krecord2t = krecord2;
				}
			}
			else
			{
				krecord2 = -1.0;
			}
			if (N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if((j + k) < krecord1)
						{
							if(pop1[curpop1][i][j][k] > 0)
							{
								krecord1 = j + k;
							}
						}
					}
				}
				if (krecord1t > krecord1)
				{
					krecord1t = krecord1;
				}
			}
			else
			{
				krecord1 = -1.0;
			}
			sprintf(str,"%s,%.2f,%.2f,%.2f,%.2f",str,N1[i],N2[i],krecord1,krecord2);
		}
	}
	else // histogram
	{
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++) // pop2
			{
				for (k=0; k<=kmax; k++)
				{
					p2hist[(j+k)] += pop2[curpop2][i][j][k];
				}
			}
			for (j=0; j<=kmax; j++) // pop1
			{
				for (k=0; k<=kmax; k++)
				{
					p1hist[(j+k)] += pop1[curpop1][i][j][k];
				}
			}
		}
	}

	if(timestep == 1)
	{
		if(krecord1t == 0)
		{
			krecord1t = -1.0;
		}
		if(krecord2t == 0)
		{
			krecord2t = -1.0;
		}
		if (krecord != 2)
		{
			fprintf(*fPointer,"%d,%d,%.2f,%.2f,%.2f,%.2f%s\n",rep+1,gen+1,N1[0],N2[0],krecord1t, krecord2t, str);
		}
		else
		{
			for (j=0; j<=2*kmax; j++)
			{
				sprintf(str,"%s,%.0f,%.0f",str,p1hist[j],p2hist[j]);
			}
			fprintf(*fPointer,"%d,%d,%.2f,%.2f%s\n",rep+1,gen+1,N1[0],N2[0], str);
		}
	}
	else
	{
		if (krecord != 2)
		{
			fprintf(*fPointer,"%d,%.2f,%.2f,%.2f,%.2f,%s\n",rep+1,N1[0],N2[0],krecord1t, krecord2t, str);
		}
		else
		{
			for (j=0; j<=2*kmax; j++)
			{
				sprintf(str,"%s,%.0f,%.0f",str,p1hist[j],p2hist[j]);
			}
			fprintf(*fPointer,"%d,%.2f,%.2f%s\n",rep+1,N1[0],N2[0], str);

		}
	}
	//sprintf(str,"");
	//free(str); 
}

// random number generating functions ran1=uniform [0,1], bnldev= binomial. gammln is needed for bnldev.
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)	
float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#define PI 3.141592654

float poidev(float xm,long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}

#undef PI

double fact(int num)
{
	// factorial
	double val = 1;
	int i;
	for (i=1; i<=num; i++)
	{
		val *= i;
	}
	return val;
}

double poipmf(double l, int k)
{
	//pmf function of poisson
	double val = (pow(l,k)*exp(-1*l))/fact(k);
	return val;
}

float gasdev(idum)
long *idum;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	float ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0 || r == 0.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}