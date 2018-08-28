/*
 Influenza Competition model 2.3
 competition model btw 1seg and 2seg ad 3SEGMENTS in a WF model
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
	int klist[8];
	int k;
};

//https://stackoverflow.com/questions/2620146/how-do-i-return-multiple-values-from-a-function-in-c
float ran1(long *seed);
float gammln(float xx);
float bnldev(float pp, int n, long *idum);
void mutate(int maxseg, long *seed, int back, int N0, double mu, int L, struct virus popop[]);
struct virus *step(int maxseg, long *seed, int rep, int t, double cost, int N0, int L, int timestep, int krecord, double s, int K, double mu, double r[],struct virus popop[],struct virus *next_gen_p,FILE **fPointer);
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int size,int a[]);

int main(int argc, char *argv[]) {
	// set progress bar and initiate timer
	clock_t begin = clock();

	// call in parameters


	char *destination = argv[1];
	char *back_s = argv[2];
	char *timestep_s = argv[3];
	char *krecord_s = argv[4];
	char *rep_s = argv[5];
	char *L_s = argv[6];
	char *s_s = argv[7];
	char *N0_s = argv[8];
	char *K_s = argv[9];
	char *mu_s = argv[10];
	char *gen_num_s = argv[11];
	char *cost_s = argv[12];
	char *N1r_s = argv[13];
	char *N2r_s = argv[14];
	char *N3r_s = argv[15];
	char *N4r_s = argv[16];
	char *N5r_s = argv[17];
	char *N6r_s = argv[18];
	char *N7r_s = argv[19];
	char *N8r_s = argv[20];
	char *seed_s = argv[21];
	char *maxseg_s = argv[22];
	char *end1;

	int back = (int) strtol(back_s,&end1,10);
	int timestep = (int) strtol(timestep_s,&end1,10);
	int krecord = (int) strtol(krecord_s, &end1, 10);
	int rep = (int) strtol(rep_s,&end1,10);
	int L = (int) strtol(L_s,&end1,10);
	double s = (double) strtof(s_s, NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	int K = (int) strtol(K_s,&end1,10);
	double mu = (double) strtof(mu_s,NULL);
	int gen_num = (int) strtol(gen_num_s,&end1,10);
	double cost = (double) strtof(cost_s,NULL);
	double Nratio[8];
	double Nratio[0] = (double) strtof(N1r_s,NULL);
	double Nratio[1] = (double) strtof(N2r_s,NULL);
	double Nratio[2] = (double) strtof(N3r_s,NULL);
	double Nratio[3] = (double) strtof(N4r_s,NULL);
	double Nratio[4] = (double) strtof(N5r_s,NULL);
	double Nratio[5] = (double) strtof(N6r_s,NULL);
	double Nratio[6] = (double) strtof(N7r_s,NULL);
	double Nratio[7] = (double) strtof(N8r_s,NULL);
	long seed = strtol(seed_s,&end1,10);
	int maxseg = (int) strtol(maxseg_s,&end1,10);	
	int i;
	double r[maxseg];
	for (i=0; i<maxseg; i++)
	{
		r[i] = (double) (pow(2,(i+1))-2)/pow(2,(i+1));
	}
	printf("back=%d, timestep=%d, krecord=%d, rep=%d, L=%d, s=%.2f, N0=%d, K=%d, mu=%.5f, gen_num=%d, cost=%.2f, N1r=%.2f, N2r=%.2f, N3r=%.2f, N4r=%.2f, N5r=%.2f, N6r=%.2f, N7r=%.2f, N8r=%.2f\n", back, timestep, krecord, rep, L, s, N0, K, mu, gen_num, cost, Nratio[0], Nratio[1], Nratio[2], Nratio[3], Nratio[4], Nratio[5], Nratio[6], Nratio[7]);
	printf("r = [%d,%d,%d,%d,%d,%d,%d,%d]\n",r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7]);

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
	sprintf(filename,"%s/c3.3s_%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f](0).csv",dest2,back,rep,L,s,N0,K,mu,gen_num,cost,Nratio[0], Nratio[1], Nratio[2], Nratio[3], Nratio[4], Nratio[5], Nratio[6], Nratio[7]);
	int filenum  = 0;
	while( access( filename, F_OK ) != -1 ) { // check if file exists and change the file number if it exists
	    filenum += 1;
		sprintf(filename,"%s/c3.3s_%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f](%d).csv",dest2,back,rep,L,s,N0,K,mu,gen_num,cost,Nratio[0], Nratio[1], Nratio[2], Nratio[3], Nratio[4], Nratio[5], Nratio[6], Nratio[7],filenum);
	}
	FILE * fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);

	// simulation start
	struct virus pop[N0];
	int N[8]; // array of subpopulation size with different segment number.
	for (i=0; i<8;i++)
	{
		N[i] = N0*Nratio[i];
	}
	struct virus next_gen[N0];
	struct virus *pop2;
	char *header = (char*) malloc(500*sizeof(char));
	if (timestep && krecord != 1) 
	{

		sprintf(header, "rep,t");
		for (i=0; i<meaxseg; i++)
		{
			sprintf(header,"%s,pop%d",header,i);		
		}
		for (i=0; i<meaxseg; i++)
		{
			sprintf(header,"%s,k%d",header,i);
		}		
		sprintf(header,"%s\n",header)
		fprintf(fPointer,header);
		free(header);
	} 
	else if (timestep && krecord == 1)
	{
		char *header1 = (char*) malloc(1250*sizeof(char));
		sprintf(header1, "rep,t,pop1,pop2,pop3");
		int m;
		for( m=0; m<51; m++) 
		{
			sprintf(header1,"%s,k1.%d",header1,m);			
		}
		for( m=0; m<51; m++)
		{
			sprintf(header1,"%s,k2.%d",header1,m);			
		}
		for( m=0; m<51; m++)
		{
			sprintf(header1,"%s,k3.%d",header1,m);
		}
		sprintf(header1,"%s\n",header1);
		fprintf("%s",header1);
		free(header1);
	}
	else
	{
		sprintf(header, "rep");
		for (i=0; i<meaxseg; i++)
		{
			sprintf(header,"%s,pop%d",header,i);		
		}
		for (i=0; i<meaxseg; i++)
		{
			sprintf(header,"%s,k%d",header,i);
		}		
		sprintf(header,"%s\n",header)
		fprintf(fPointer,header);
		free(header);
	}
	int j,add,repe,gen;
	for (repe=0; repe<rep; repe++)
	{	
		// fill in the pop array
		add = 0; // necessary for indexing pop.
		for (j=0; j<mexseg; j++)
		{
			if (j!= 0)
			{
				add += N[j-1];		
			}
			for (i=0;i<N[j];i++)
			{
				pop[i + add].id = j;
			}
		}
		printf("REP=%d/%d\n",repe,rep);
		for (gen=0; gen<gen_num; gen++)
		{
			// run through generation
			//printf("GEN=%d/%d\n",gen,gen_num);
			mutate(maxseg,&seed,back,N0,mu,L,pop);
			pop2 = step(maxseg,&seed,(repe+1),(gen+1),cost,N0,L,timestep,krecord,s,K,mu,r,pop,next_gen,&fPointer);
			memcpy(pop,pop2,sizeof(struct virus)*N0); // cycle between pop and pop2 to continue looping.
		}
	}

	// close file and timer
	fclose(fPointer);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent was %.2f minutes\n", time_spent/60.0);

	return 0;
}


struct virus *step(int maxseg,long *seed,int rep, int t, double cost, int N0, int L, int timestep, int krecord, double s, int K, double mu, double r[], struct virus popop[],struct virus *next_gen_p,FILE **fPointer) {
//struct virus *step(int rep, int t,struct virus popop[],struct virus *next_gen_p,FILE **fPointer) {
	// goes through reproduction process
	// the process is depicted at the top of the script.
	// input: pop struct array
	// output: next generation's pop struct array
	int l = 0; // next gen length
	int ks[maxseg];
	int ksl[maxseg];
	int s1, s2;
	int i,j;
	int mink[maxseg]; // array of minimum k value for different segmented viruses
	int ksum; // when 2 parents have same # of segments, add all the k's from random combinations to determine if the baby will survive
	for (i=0; i<maxseg; i++)
	{
		mink[i] = 1000;
	}
	float random;
	while (l < N0) 
	{
		s1 = floor(ran1(seed)*N0); // sample 1
		s2 = floor(ran1(seed)*N0); // sample 2
		//printf("s1:%d k1:%d k2:%d k:%d\n",s1,popop[s1].k1,popop[s1].k2,popop[s1].k);
		//printf("s2:%d k1:%d k2:%d k:%d\n",s2,popop[s2].k1,popop[s2].k2,popop[s2].k);
		if ((popop[s1].id == 1 && popop[s2].id == 1) || (popop[s1].id != popop[s2].id)) 
		{
			// either parents is segment 1 or both parents are differently segmented
			if (ran1(seed) < 0.5) 
			{	
				// pick s1
				if (ran1(seed) < pow(1.0-s,popop[s1].k)) 
				{
					next_gen_p[l].id = popop[s1].id;
					for (i=0; i<next_gen_p[l].id;i++)
					{
						next_gen_p[l].klist[i] = popop[s1].klist[i];
					}
					next_gen_p[l].k = popop[s1].k;
					if (krecord == 0)
					{
						ks[next_gen_p[l].id] += next_gen_p[l].k
						ksl[next_gen_p[l].id]++;
					}
					else if (krecord == 2)
					{
						if (ks[next_gen_p[l].id] < mink[next_gen_p[l].id])
						{
							ks[next_gen_p[l].id] = next_gen_p[l].k
							ksl[next_gen_p[l].id]++;							
						}
					}
					l++;
				}
			} 
			else 
			{ 
				// pick s2
				if (ran1(seed) < pow(1.0-s,popop[s2].k))
				{
					next_gen_p[l].id = popop[s1].id;
					for (i=0; i<next_gen_p[l].id;i++)
					{
						next_gen_p[l].klist[i] = popop[s1].klist[i];
					}
					next_gen_p[l].k = popop[s1].k;
					if (krecord == 0)
					{
						ks[next_gen_p[l].id] += next_gen_p[l].k
						ksl[next_gen_p[l].id]++;
					}
					else if (krecord == 2)
					{
						if (ks[next_gen_p[l].id] < mink[next_gen_p[l].id])
						{
							ks[next_gen_p[l].id] = next_gen_p[l].k
							ksl[next_gen_p[l].id]++;							
						}
					}
					l++;
				}
			}
		}
		else
		{ // both parents are same segmented
			next_gen_p[l].id = popop[s1].id;
			for (j=0; j<next_gen_p[l].id;j++)
			{
				if (ran1(seed) < 0.5)
				{
					next_gen_p[l].klist[j] = popop[s1].klist[j];
					ksum += popop[s1].klist[j]
				}
				else
				{
					next_gen_p[l].klist[j] = popop[s2].klist[j];
					ksum += popop[s2].klist[j]
				}
			}
			if (ran1(seed) < pow(1.0-s,ksum)*(1.0-cost))
			{
				if (krecord == 0)
				{
					ks[next_gen_p[l].id] += next_gen_p[l].k
					ksl[next_gen_p[l].id]++;
				}
				else if (krecord == 2)
				{
					if (ks[next_gen_p[l].id] < mink[next_gen_p[l].id])
					{
						ks[next_gen_p[l].id] = next_gen_p[l].k
						ksl[next_gen_p[l].id]++;							
					}
				}
				ksum = 0;
				l++;
			}
			else
			{
				ksum = 0;
			}
		}			
	}
	char *line = (*char) malloc(300*sizeof(char));
	char *line1 = (*char) malloc(150*sizeof(char));
	char *line2 = (*char) malloc(150*sizeof(char));
	if (krecord == 0)
	{
		float mk[maxseg];
		for (i=0; i<maxseg; i++)
		{
			if (ksl[i] == 0)
			{
				mk[i] = -1;
			}
			else
			{
				mk[i] = (float) ks[i]/ksl[i];
			}
			sprintf("%s,%.2f",line1,mk[i]);
			sprintf("%s,%d",line2,ksl[i]);
		}
		sprintf(line,"%d,%d%s%s\n",rep,t,line1,line2);
		fprintf(*fPointer,line);
	} 
	else if (krecord == 1) 
	{
	/*	int k1list[51];
		int k2list[51];
		int k3list[51];
		for (i=0; i<ks1l; i++) 
		{
			if (ks1[i] > 50) 
			{
				printf("WARNINGWARNING THERE'S INDIVIDUAL WITH K HIGHER THAN 50!!!!\n");
			}
			k1list[ks1[i]] += 1; 
		}
		for (i=0; i<ks2l; i++)
		{
			if (ks2[i] > 50) 
			{
				printf("WARNINGWARNING THERE'S INDIVIDUAL WITH K HIGHER THAN 50!!!!\n");
			}
			k2list[ks2[i]] += 1;
		}
		for (i=0; i<ks3l; i++)
		{
			if (ks3[i] > 50) 
			{
				printf("WARNINGWARNING THERE'S INDIVIDUAL WITH K HIGHER THAN 50!!!!\n");
			}
			k3list[ks3[i]] += 1;
		}
		char *line = (char*) malloc(20100*sizeof(char));
		sprintf(line,"%d,%d,%d,%d,%d\n",rep,t,ks1l,ks2l,ks3l);
		for (i=0; i<51; i++)
		{
			sprintf(line,"%s,%d",line,k1list[i]);
		}
		for (i=0; i<51; i++)
		{
			sprintf(line,"%s,%d",line,k2list[i]);
		}
		for (i=0; i<51; i++)
		{
			sprintf(line,"%s,%d",line,k3list[i]);
		}
		fprintf(*fPointer,"%s\n",line);
		free(line);
	*/
	}
	else
	{
		for (i=0; i<maxseg; i++)
		{
			if (ksl[i] == 0)
			{
				ks[i] = -1;
			}
			sprintf("%s,%d",line1,ks[i]);
			sprintf("%s,%d",line2,ksl[i]);
		}
		sprintf(line,"%d,%d%s%s\n",rep,t,line1,line2);
		fprintf(*fPointer,line);
	}	
	free(line);
	free(line1);
	free(line2);
	return next_gen_p;
}


void mutate(int maxseg, long *seed, int back, int N0, double mu, int L, struct virus popop[]) {
	// goes through population and make every inidividual go through 
	// mutation process.
	// input: population struct array.
	// output: population struct with updated k values according to
	//   mutation rate.
	int i,j,k,mut_num,ksum;
	float random;
	int random2;
	for (i=0;i<N0; i++)
	{
		mut_num = bnldev(mu,L,seed); // binomial pick of number of mutation based on mu.
		for (j=0; j<mut_num; j++)
		{
			if (back == 1)
			{
				// back mutation
				if (ran1(seed)*L < popop[i].k)
				{
					// back mutation happens
					random = ran1(seed);
					for (k=0; k<popop[i].id; k++)
					{
						ksum += popop[i].klist[k];
							if (random < (float) ksum/popop[i].k)
							{
								popop[i].klist[k]--;
								popop[i].k--;
							}
					}
				}
				else
				{
					// back mutation doesn't happen
					random2 = floor(ran1(seed)*maxseg);
					popop[i].klist[random2]++;
					popop[i].k++;
				}
			}
			else
			{
				// no back mutation
				random2 = floor(ran1(seed)*maxseg);
				popop[i].klist[random2]++;
				popop[i].k++;
			}
		}
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

int intmin(int argc, int array[]){
  int min_num = array[0];
  int i;
  for(i=1;i<argc;i++)
  {
    if(array[i]<min_num)
    {
      min_num = array[i];
    }
  }
  return min_num;
}

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++)
  {
    sum += a[i];
  }
  return sum;
}
