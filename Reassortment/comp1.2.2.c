/*
Influenza Competition model 1.2.2
Same concept as comp1.2, but instead of deciding number of progenies through poisson process with rates related to fiteness(w),
it randomly picks 2 parents and decides if that progeny will survive with probability related to fitness.

Some key differences:
- back mutation, tracking info, and program cut strategy parameterized.
- No sequence information using array. Only keeps track of mutant allele amount
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

/*
 Parameters
 back = whether the simulation allows back mutation (1 if there's one. 0 if there isn't)
 timestep = whether a simulation wants to record every generation or only the end generation per repetition
 krecord = how you want to record k in the output csv. 0= mean k. 1=all k values of indivs in an array. 2= minimum k.
 N0 = Population size
 K = Carrying capacity
 L = sequence length
 s = fitness decrease from deleterious mutation
 mu = mutation rate per site
 gen_num = generation amount 
 r = reassortment rate
 rep = repetition amount
 N1r = ratio of 1segment virus
*/

// virus basic structure
struct virus {
	int id;
	int k1;
	int k2;
	int k;
	int progeny;
};

float ran1(long *seed);
float gammln(float xx);
float bnldev(float pp, int n, long *idum);
float poidev(float xm,long *idum);
void mutate(long *seed, int back, int N, double mu, int L, struct virus popop[]);
struct virus *step(long *seed,int rep, int t, double cost, int N, int L, int timestep, int krecord, double s, int K, double mu, double r,struct virus popop[],struct virus *next_gen_p,FILE **fPointer,int* N1, int* N2, int gen, int gen_num);
int *intsample(int *array, int l, int n, int replacement, long *seed, int *receive);
int *intdel(int *array, int l, int ran, int *receive);
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int size,int a[]);

int main(int argc, char *argv[]) {
	// set progress bar and initiate timer
	clock_t begin = clock();
	
	char *destination = argv[1];
	char *back_s = argv[2];
	char *timestep_s = argv[3];
	char *krecord_s = argv[4];
	char *untilext_s = argv[5];
	char *rep_s = argv[6];
	char *L_s = argv[7];
	char *s_s = argv[8];
	char *N0_s = argv[9];
	char *K_s = argv[10];
	char *mu_s = argv[11];
	char *gen_num_s = argv[12];
	char *cost_s = argv[13];
	char *r_s = argv[14];
	char *N1r_s = argv[15];
	char *seed_s = argv[16];
	char *end1;
	
	int back = (int) strtol(back_s,&end1,10);
	int timestep = (int) strtol(timestep_s,&end1,10);
	int krecord = (int) strtol(krecord_s, &end1,10);
	int untilext = (int) strtol(untilext_s,&end1,10);
	int rep = (int) strtol(rep_s,&end1,10);
	int L = (int) strtol(L_s,&end1,10);
	double s = (double) strtof(s_s,NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	int K = (int) strtol(K_s,&end1,10);
	double mu = (double) strtof(mu_s,NULL);
	int gen_num = (int) strtol(gen_num_s,&end1,10);
	double cost = (double) strtof(cost_s,NULL);
	double r = (double) strtof(r_s,NULL);
	double N1r = (double) strtof(N1r_s,NULL);
	long seed = strtol(seed_s,&end1,10);

	/*
	int back = 0;
	int timestep = 1;
	int krecord = 0;
	int untilext = 0;
	int rep = 1;
	int L = 300;
	double s = 0.05;
	int N0 = 10000;
	int K = 10000;
	double mu = 0.06;
	int gen_num = 2;
	double cost = 0.0;
	double r = 0.5;
	double N1r = 0.5;
	char *destination = "ctest";
	long seed = -2389;
	*/
	printf("back=%d, timestep=%d, krecord=%d, untilext=%d, rep=%d, L=%d, s=%.2f, N0=%d, K=%d, mu=%.5f, gen_num=%d, cost=%.2f, r=%.2f, N1r=%.2f\n", back, timestep, krecord, untilext, rep, L, s, N0, K, mu, gen_num, cost, r, N1r);
		
	//initiate csv file
	//// set up folder
	char *dest2 = (char*) malloc(50*sizeof(char)); 
	sprintf(dest2, "./data/%s", destination);
	struct stat st = {0};
	if (stat(dest2, &st) == -1) { // if the destination folder doesn't exist, make one with that name.
    	mkdir(dest2, 0700);
	}
	//// initiate file
	char *filename = (char*) malloc(100*sizeof(char));
	sprintf(filename,"%s/c1.2.2s_%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f(0).csv",dest2,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r);
	int filenum  = 0;
	while( access( filename, F_OK ) != -1 ) { // check if file exists and change the file number if it exists
	    filenum += 1;
		sprintf(filename,"%s/c1.2.2s_%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f(%d).csv",dest2,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r,filenum);
	}
	FILE * fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);

	// simulation start
	if (timestep && krecord != 1) 
	{
		fprintf(fPointer,"rep,t,pop1,pop2,k1,k2\n");
	} 
	else if (timestep && krecord == 1)
	{
		char *header1 = (char*) malloc(650*sizeof(char));
		sprintf(header1, "rep,t,pop1,pop2,");
		int m;
		for( m=0; m<51; m++) 
		{
			sprintf(header1,"%s,k1.%d",header1,m);			
		}
		for( m=0; m<51; m++)
		{
			sprintf(header1,"%s,k2.%d",header1,m);			
		}
		sprintf(header1,"%s\n",header1);
		printf("%s",header1);
		free(header1);
	}
	else
	{
		fprintf(fPointer,"pop1,pop2,k1,k2\n");
	}

	// initiate some variables that will change over generations
	int N1; // seg1 population
	int N2; // seg2 population
	int N; // N1 + N2
	struct virus pop1[K+K/2];
	struct virus next_gen[K+K/2];
	struct virus *pop2;
	int repe,i,gen;
	for (repe=0; repe<rep; repe++){
		//if (repe % 100 == 0)
		{
			printf("REP=%d/%d\n",repe,rep);
		}
		// initialize pop (generation 0)
		N1 = N0*N1r; // initial 1seg pop
		N2 = N0*(1-N1r); // initial 2seg pop
		N = N0;
		for (i=0;i<N1;i++){
			pop1[i].id = 1;
			pop1[i].k1 = 0;
			pop1[i].k2 = 0;
			pop1[i].k = 0;
			pop1[i].progeny = 0;
		}
		for (i=0;i<N2;i++){
			pop1[i+N1].id = 2;
			pop1[i+N1].k1 = 0;
			pop1[i+N1].k2 = 0;
			pop1[i+N1].k = 0;		
			pop1[i+N1].progeny = 0;

		}
		for (gen=0; gen<gen_num; gen++){ // run through generation
			//printf("GEN=%d/%d\n",gen+1,gen_num);
			// cycle btw pop and popb to continue looping.
			mutate(&seed, back, N, mu, L, pop1); // seg1's mutate
			pop2 = step(&seed,(repe+1),(gen+1),cost,N,L,timestep,krecord,s,K,mu,r,pop1,next_gen,&fPointer,&N1,&N2,gen,gen_num);
			//printf("worked til gen=%d\n",gen);
			//printf("pop1[0].progeny = %d\n",pop1[0].progeny);
			//printf("pop2[0].progeny = %d\n",pop2[0].progeny);
			memcpy(pop1,pop2,sizeof(struct virus)*N); // cycle between pop and pop2 to continue looping.			
			//printf("pop1[0].progeny = %d\n",pop1[0].progeny);
			//printf("pop2[0].progeny = %d\n",pop2[0].progeny);

			N = N1 + N2;
		}
	}
	// close file and timer
	fclose(fPointer);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent was %.2f minutes\n", time_spent/60.0);
	return 0;
}

struct virus *step(long *seed,int rep, int t, double cost, int N, int L, int timestep, int krecord, double s, int K, double mu, double r,struct virus popop[],struct virus *next_gen_p,FILE **fPointer,int* N1, int* N2, int gen, int gen_num)
{
	// goes through reproduction process
	// the process is depicted at the top of the script.
	// input: pop struct array
	// output: next generation's pop struct array
	int remaining[K+K/2]; // list of viruses who still have progenise to give.
	int remainingl = 0;
	int i;
	int l = 0; // length of next_gen_p


	for (i=0;i<N;i++) // decide how many youngs each virus will have and fill the 'remaining' array.
	{
		popop[i].progeny = poidev((2/(1+(float)N/K)),seed);
		if (popop[i].progeny > 0)
		{
			remaining[remainingl] = i;
			remainingl++;
		}
	}

	//printf("remainings: ");
	//for (i=0; i<remainingl; i++) {
	//	printf("%d, ",remaining[i]);
	//}
	int *samp;
	int sampreceive[2];
	int temp;
	int *temparray;
	int remainreceive[remainingl];
	int sampindex; //samp index in remaining;
	float kval1;
	int kval1l = 0;
	float kval2;
	int kval2l = 0;
	int sampturns; // used for for-loop variable
	int sampturns2;
	if (krecord == 2)
	{
		kval1 = 100000;
		kval2 = 100000;
	}
	while (remainingl >= 2)
	{
		samp = intsample(remaining,remainingl, 2, 0, seed, sampreceive); // sample 2 parent viruses
		if (samp[0] < samp[1]) {
			temp = samp[0];
			samp[0] = samp[1];
			samp[1] = temp;
		}
		//printf("sample[0]=%d\n",samp[0]);	
		//printf("sample[1]=%d\n",samp[1]);
		if ( (popop[samp[0]].id == 1 || popop[samp[1]].id == 1) || (ran1(seed) < r) )
		{
			for (sampturns=0; sampturns < 2; sampturns++)
			{
				if (ran1(seed) < pow(1.0-s,popop[samp[sampturns]].k))
				{					
					next_gen_p[l].id = popop[samp[sampturns]].id;
					next_gen_p[l].k1 = popop[samp[sampturns]].k1;
					next_gen_p[l].k2 = popop[samp[sampturns]].k2;
					next_gen_p[l].k = popop[samp[sampturns]].k;
					if (krecord == 0)
					{
						if (next_gen_p[l].id == 1)
						{
							kval1 += next_gen_p[l].k; 
							kval1l++;
						}
						else
						{
							kval2 += next_gen_p[l].k;
							kval2l++;
						}
					}
					else if (krecord == 2)
					{
						if (next_gen_p[l].id == 1 && next_gen_p[l].k < kval1)
						{
							kval1 = next_gen_p[l].k; 
							kval1l++;
						}
						else if (next_gen_p[l].id == 2 && next_gen_p[l].k < kval2)
						{
							kval2 += next_gen_p[l].k;
							kval2l++;
						}
					}
					l++;
				}
				popop[samp[sampturns]].progeny -= 1;
			}	
		}
		else
		{
			for (sampturns=0; sampturns < 2; sampturns++)
			{
				if (sampturns == 0)
				{
					sampturns2 = 1;
				}
				else
				{
					sampturns2 = 0;
				}
				if (ran1(seed) < pow(1.0-s,popop[samp[sampturns]].k1 + popop[samp[sampturns2]].k2))
				{					
					next_gen_p[l].id = popop[samp[sampturns]].id;
					next_gen_p[l].k1 = popop[samp[sampturns]].k1;
					next_gen_p[l].k2 = popop[samp[sampturns2]].k2;
					next_gen_p[l].k = next_gen_p[l].k1 + next_gen_p[l].k2;
					if (krecord == 0)
					{
							kval2 += next_gen_p[l].k;
							kval2l++;
					}
					else if (krecord == 2)
					{
						if (next_gen_p[l].k < kval2)
						{
							kval2 += next_gen_p[l].k;
							kval2l++;
						}
					}
					l++;
				}	
				popop[samp[sampturns]].progeny -= 1;
			}	
		}
		for (sampturns=0; sampturns < 2; sampturns++)
		{
			if (popop[samp[sampturns]].progeny == 0)
			{
				for (i= (remainingl-1); i>-1; i--)
				{
					if (samp[sampturns] == remaining[i])
					{
						sampindex = i;
					}

				}
				temparray = intdel(remaining,remainingl,sampindex,remainreceive);
				//printf("sample[0] %d is erased!\n",samp[0]);
				remainingl--;
				memcpy(remaining, temparray, sizeof(int)*remainingl);
			}			
		}

		/*
		printf("remainings: ");
		for (i=0; i<remainingl; i++) {
			printf("%d, ",remaining[i]);
		}
		printf("\n");
		*/
	}
	if (l - 1 > K+K/2) 
	{
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNING: NEXT_GEN ALLOCATION WAS NOT ENOUGH. ALLOCATE MORE MEMORIES TO RUN THIS SIMULATION\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
	}
	if (remainingl > 0)
	{				
		for (i=0; i<popop[remaining[0]].progeny; i++) 
		{
			if (ran1(seed) < pow(1.0-s,popop[remaining[0]].k))
			{
				next_gen_p[l].id = popop[remaining[0]].id;
				next_gen_p[l].k1 = popop[remaining[0]].k1;
				next_gen_p[l].k2 = popop[remaining[0]].k2;
				next_gen_p[l].k = popop[remaining[0]].k;
			}
			if (krecord == 0)
			{
				if (next_gen_p[l].id == 1)
				{
					kval1 += next_gen_p[l].k;
					kval1l++; 
				}
				else
				{
					kval2 += next_gen_p[l].k;
					kval2l++;
				}
			}
			else if (krecord == 2)
			{
				if (next_gen_p[l].id == 1 && next_gen_p[l].k < kval1)
				{
					kval1 = next_gen_p[l].k; 
					kval1l++;
				}
				else if (next_gen_p[l].id == 2 && next_gen_p[l].k < kval2)
				{
					kval2 += next_gen_p[l].k;
					kval2l++;
				}
			}
			l++;					
		}
	}

	/*for (i=0; i<*N2; i++) {
		printf("pop%d has k=%d\n",i,next_gen_p[i].k);
	}*/
	*N1 = kval1l;
	*N2 = kval2l;
	if (krecord == 0)
	{
		if (kval1l == 0) 
		{
			kval1 = -1.0;
		} 
		else 
		{
			kval1 = (float) kval1/kval1l;
		}
		if (kval2l == 0) 
		{
			kval2 = -1.0;
		} 
		else 
		{
			kval2 = (float) kval2/kval2l;
		}
	}
	else if (krecord == 1) 
	{
		// WILL FINISH WHEN NECESSARY.
		/*if (ks1l == 0){
			char *k1str = "'NA'"
		} else {
			char *k1str	= (char *) malloc(N0*4*sizeof(char));
		}
		if (ks2l == 0){
			char *k2str = "'NA'"
		} else {
			char *k2str = (char *) malloc(N0*4*sizeof(char));
			ks2strl = 0;
			for (int i=0; i<ks2l; i++) {
				char *int2str = (char *) malloc(4);
				sprintf(int2str,"%d",ks2[i]);
				k2str[i] 
			}
		}
		fprintf(fPointer,"",);*/
	}
	else 
	{
		if (kval1l == 0) 
		{
			kval1 = -1.0;
		} 

		if (kval2l == 0)
		{
			kval2 = -1.0;
		}
	}
	if (timestep)
		{
			fprintf(*fPointer,"%d,%d,%d,%d,%.2f,%.2f\n",rep,t,kval1l,kval2l,kval1,kval2);
		}
		else if (gen == (gen_num -1) || kval1l == 0 || kval2l == 0)
		{
			fprintf(*fPointer,"%d,%d,%.2f,%.2f\n",kval1l,kval2l,kval1,kval2);
		}
	//printf("-----worked til gen#=%d--------\n",gen);
	return next_gen_p;	
}


int *intsample(int *array, int l, int n, int replacement,long *seed, int *receive) {
	// sampling from an integer array
	// array = list of integer to sapmle from
	// l = list length
	// n = sample size (don't put more than the length of list if no replacement)
	// replacement = true or false for replacement
	int i, ll, ran;
	ll = l;
	if (replacement)
	{
		for( i=0; i<n; i++)
		{
			ran = floor(ran1(seed)*ll);
			receive[i] = array[ran];
		}
	}
	else
	{
		int receive2[ll];
		for( i=0; i<n; i++)
		{
			ran = floor(ran1(seed)*ll);
			receive[i] = array[ran];
			array = intdel(array,ll,ran,receive2);
			ll -= 1;
		}
	}
	return receive;
}

int *intdel(int *array, int l, int ran, int *receive)
{
	// deletes an specified element in an array and returns the revised one.
	// array = an array where an element will be deleted
	// l = array length
	// ran = element index in an array

	int i;
	for (i=0; i<ran; i++)
	{
		receive[i] = array[i];
	}
	for (i=ran+1; i<l; i++)
	{
		receive[i-1] = array[i];
	}
	return receive;
}


void mutate(long *seed, int back, int N, double mu, int L, struct virus popop[]) {
	// goes through population and make every inidividual go through 
	// mutation process.
	// input: population struct array.
	// output: population struct with updated k values according to
	//   mutation rate.
	int i,j,mut_num,breakpt;
	for (i=0;i<N; i++)
	{
		mut_num = bnldev(mu,L,seed); // binomial pick of number of mutation based on mu.
		if (back == 1)
		{ // back mutation
			if (popop[i].id == 1) 
			{ // individual is 1segment
				for (j=0; j<mut_num; j++) 
				{ // add or subtract k
					if (ran1(seed)*L < popop[i].k) 
					{ // k/L is the probability of back mutating.
						popop[i].k -= 1;
					}
					else
					{
						popop[i].k += 1;
					}
				}
			} 
			else
			{ // individual is 2segment
				for (j=0; j<mut_num; j++) 
				{ // add or subtract k1 or k2, and k.
					if (ran1(seed)*L < popop[i].k) 
					{
						if (popop[i].k1==0) 
						{
							popop[i].k2 -= 1;
							popop[i].k -= 1;
						} 
						else if (popop[i].k2 == 0) 
						{
							popop[i].k1 -= 1;
							popop[i].k -= 1;
						} 
						else 
						{
							if (ran1(seed) < (float)popop[i].k1/popop[i].k) 
							{
								popop[i].k1 -= 1;
								popop[i].k -= 1;
							} 
							else 
							{
								popop[i].k2 -= 1;
								popop[i].k -= 1;
							}
						}
					}
					else 
					{
						if (ran1(seed) < 0.5) 
						{
							popop[i].k1 += 1;
							popop[i].k += 1;
						}
						else
						{
							popop[i].k2 += 1;
							popop[i].k += 1;
						}
					}
				}
			}
		}
		else 
		{ // no back mutation 
			if (popop[i].id == 1) 
			{
				popop[i].k += mut_num;
			}
			else 
			{
				breakpt = floor(ran1(seed)*(mut_num+1));
				popop[i].k1 += breakpt;
				popop[i].k2 += mut_num - breakpt;
				popop[i].k += mut_num;
			}
		}
		//printf("%d,",popop[i].k);
	}
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
float bnldev(float pp, int n,long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) bnl += 1.0;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI

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

int intmin(int argc, int array[]){
  int min_num = array[0];
  int i;
  for(i=1;i<argc;i++){
    if(array[i]<min_num){
      min_num = array[i];
    }
  }
  return min_num;
}

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++){
    sum += a[i];
  }
  return sum;
}
