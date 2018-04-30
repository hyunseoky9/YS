/*
 Influenza Competition model 1.3
 competition model btw 1seg and 2seg in a WF model
 Some key differences from 1.2's:
 - WF model (no demographic factor or carrying capcacity for that matter.)
 Reproduction process:
 1. pick 2 random parents
 2. If either one is 1segment, only one of them gets replicated.
 3. If both of them are 2segement, it goes through recombination with 
    a probability of r, and one of the progeny is created.
 4. the created progeny from either proces  2. or 3. are assessed 
 5. whether they will survive by their fitness factor w.
 6. If a random number ([0,1]) exceeds w, the progeny dies and 
    has to go through the process 2.-5. again.
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
	int k1;
	int k2;
	int k;
};

//https://stackoverflow.com/questions/2620146/how-do-i-return-multiple-values-from-a-function-in-c
float ran1(long *seed);
float gammln(float xx);
float poidev(float xm,long *idum);
void mutate(long *seed, int krecord, int N1orN2, int N0, float initkmu, float *beforek, struct virus pop[]);
void step(long *seed, int N1orN2, int N0, int krecord, double s, double r,struct virus pop[], float *afterk);
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int size,int a[]);

int main(int argc, char *argv[]) {
	// set progress bar and initiate timer
	clock_t begin = clock();

	// call in parameters
	
	char *destination = argv[1];
	char *krecord_s = argv[2];
	char *rep_s = argv[3];
	char *s_s = argv[4];
	char *N0_s = argv[5];
	char *initkmu_s = argv[6];
	char *cost_s = argv[7];
	char *r_s = argv[8];
	char *N1orN2_s = argv[9]; // 1 if only 1segments; 2 if only 2 segments
	char *seed_s = argv[10];
	char *end1;

	int krecord = (int) strtol(krecord_s, &end1, 10);
	int rep = (int) strtol(rep_s,&end1,10);
	double s = (double) strtof(s_s, NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	float initkmu = strtof(initkmu_s,NULL);
	double cost = (double) strtof(cost_s,NULL);
	double r = (double) strtof(r_s,NULL);
	int N1orN2 = (int) strtol(N1orN2_s,&end1,10);
	long seed = strtol(seed_s,&end1,10);

	printf("destination=%s, krecord=%d, rep=%d, s=%.2f, N0=%d, initkmu=%.5f, cost=%.2f, r=%.2f, N1orN2=%d\n", destination, krecord, rep, s, N0, initkmu , cost, r, N1orN2);
	
	/*
	char *destination = "ctest";
	int krecord = 0;
	int rep = 1;
	double s = 0.05;
	int N0 = 12;
	float initkmu = 2;
	double cost = 0.0;
	double r = 0.5;
	int N1orN2 = 1;
	long seed = -235;
	*/


	//initiate csv filename
	//// set up folder
	char *dest2 = (char*) malloc(50*sizeof(char)); 
	sprintf(dest2, "./data/%s", destination);
	struct stat st = {0};
	if (stat(dest2, &st) == -1) { // if the destination folder doesn't exist, make one with that name.
    	mkdir(dest2, 0700);
	}
	//// initiate file
	char *filename = (char*) malloc(100*sizeof(char));
	sprintf(filename,"%s/step_%d,%.2f,%d,%.3f,%.2f,%.2f,%d(0).csv",dest2,rep,s,N0,initkmu,cost,r,N1orN2);
	int filenum  = 0;
	while( access( filename, F_OK ) != -1 ) { // check if file exists and change the file number if it exists
	    filenum += 1;
		sprintf(filename,"%s/step_%d,%.2f,%d,%.3f,%.2f,%.2f,%d(%d).csv",dest2,rep,s,N0,initkmu,cost,r,N1orN2,filenum);
	}	
	FILE * fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);

	// simulation start
	struct virus pop[N0];
	fprintf(fPointer,"rep,k0,k1,k1-k0\n");
	int repe;
	float beforek, afterk;
	float kdistmean = 0;
	for (repe=0; repe<rep; repe++){	
		/*printf("REP=%d/%d\n",repe+1,rep);
		for (int j=0; j<N0; j++)
		{
			printf("pop[%d] k=%d\n",j,pop[j].k);
		}*/
		mutate(&seed, krecord, N1orN2, N0, initkmu, &beforek, pop);
		//printf("--------------after mutation -------------\n");

		step(&seed,N1orN2,N0,krecord,s,r,pop,&afterk);
		fprintf(fPointer,"%d,%f,%f,%f\n",repe+1,beforek,afterk,afterk-beforek);
		kdistmean += (afterk-beforek);
	}
	fprintf(fPointer,"%.3f,0,0,0",kdistmean);
	kdistmean = kdistmean/(float)rep;
	// close file and timer
	fclose(fPointer);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent was %.2f minutes\n", time_spent/60.0);
	printf("k distance mean is: %f\n",kdistmean);
	return 0;
}


void step(long *seed, int N1orN2, int N0, int krecord, double s, double r,struct virus pop[], float *afterk) {

	// goes through reproduction process
	// the process is depicted at the top of the script.
	// input: pop struct array
	// output: next generation's pop struct array
	int l = 0; // next gen length
	int ks[N0];
	int ksl = 0;
	int s1, s2;
	if (N1orN2 == 1)
	{
		while (l < N0) 
		{
			s1 = floor(ran1(seed)*N0); // sample 1
			if (ran1(seed) < pow(1.0-s,pop[s1].k))
			{
				ks[ksl] = pop[s1].k;
				ksl += 1;
				l += 1;
			}				
		}
	}
	else
	{
		while (l < N0) 
		{
			s1 = floor(ran1(seed)*N0); // sample 1
			s2 = floor(ran1(seed)*N0); // sample 2
			if (ran1(seed)<r) { // recombination happens
				if (ran1(seed) < 0.5) { // pick k1 from s1 and k2 from s2
					if (ran1(seed) < pow(1.0-s,(pop[s1].k1 + pop[s2].k2))){
						ks[ksl] = pop[s1].k1 + pop[s2].k2;
						ksl += 1;
						l += 1;
					}
				} else { // pick k1 from s2 and k2 from s1
					if (ran1(seed) < pow(1.0-s,(pop[s1].k2 + pop[s2].k1))){
						ks[ksl] = pop[s1].k2 + pop[s2].k1;
						ksl += 1;
						l += 1;
					}
				}
			} else { // recombintion doesn't happen
				if (ran1(seed) < 0.5) { // pick s1
					if (ran1(seed) < pow(1.0-s,pop[s1].k)){
						ks[ksl] = pop[s1].k;
						ksl += 1;
						l += 1;
					}
				} else { // pick s2
					if (ran1(seed) < pow(1.0-s,pop[s2].k)){
						ks[ksl] = pop[s2].k;
						ksl += 1;				
						l += 1;
					}
				}
			}
				
		}
	}
	float mk;
	if (krecord == 0){
		if (ksl == 0) {
			mk = -1.0;
		} else {
			mk = (float)intsum(ksl,ks)/ksl;		
		}
	} else if (krecord == 1) {
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
	} else {
		if (ksl == 0) {
			mk = -1.0;
		} else {
			mk = intmin(ksl,ks);
		} 
	}
	*afterk = mk;
	/*for(int j=0; j<N0; j++)
	{
		printf("ks[%d] k=%d\n",j,ks[j]);
	}
	printf("mk in step= %f\n",*afterk);*/
}

void mutate(long *seed, int krecord, int N1orN2, int N0, float initkmu, float *beforek, struct virus pop[]) {
	// goes through population and make every inidividual go through 
	// mutation process.
	// input: population struct array.
	// output: population struct with updated k values according to
	//   mutation rate.
	/*printf("\n\n");
		for (int j=0; j<N0; j++)
	{
		printf("pop[%d] k=%d\n",j,pop[j].k);
	}*/
	int i,mut_num,breakpt;
	int ks[N0];
	int ksl = 0;
	float mk;
	if (N1orN2 == 1){
		for (i=0;i<N0; i++){
			mut_num = poidev(initkmu,seed); // binomial pick of number of mutation based on mu.
			pop[i].k = mut_num;
			ks[ksl] = pop[i].k;
			ksl++;		
		}
	}
	else
	{
		for (i=0;i<N0; i++){
			mut_num = poidev(initkmu,seed); // binomial pick of number of mutation based on mu.
			breakpt = floor(ran1(seed)*(mut_num+1));
			pop[i].k1 = breakpt;
			pop[i].k2 = mut_num - breakpt;
			pop[i].k = mut_num;
			ks[ksl] = pop[i].k;
			ksl++;
		}
	}

	if (krecord == 0){
		if (ksl == 0) {
			mk = -1.0;
		} else {
			mk = (float)intsum(ksl,ks)/ksl;		
		}
	} else if (krecord == 1) {
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
	} else {
		if (ksl == 0) {
			mk = -1.0;
		} else {
			mk = intmin(ksl,ks);
		} 
	}
	/*
	printf("\n\n");
		for (int j=0; j<N0; j++)
	{
		printf("pop[%d] k=%d\n",j,pop[j].k);
	}
	printf("\n\n");*/
	*beforek = mk;
	//printf("mk = %f\n",*beforek);
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
