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

int back = 0;
int timestep = 0;
int krecord = 0;
char *destination = "ctest";
int rep = 1;
int L = 300;
double s = 0.05;
#define N0 50
int K = 1000;
float mu = 0.005;
int gen_num = 10;
double cost = 0.00;
double r = 0.5;
double N1r = 0.5;
long seed = -1;

struct virus {
	int id;
	int k1;
	int k2;
	int k;
};

//https://stackoverflow.com/questions/2620146/how-do-i-return-multiple-values-from-a-function-in-c
float ran1(long *seed);
float gammln(float xx);
float bnldev(float pp, int n, long *idum);
struct virus *step(struct virus popop[],struct virus *next_gen_p);


int main(void) {
	//initialize pop
	struct virus pop[N0];
	int N1 = N0*N1r;
	int N2 = N0*(1-N1r);
	for (int i=0;i<N1;i++){
		pop[i].id = 1;
		pop[i].k1 = 0;
		pop[i].k2 = 0;
		pop[i].k = 0;
	}
	for (int i=0;i<N2;i++){
		pop[i+N1].id = 2;
		pop[i+N1].k1 = 0;
		pop[i+N1].k2 = 0;
		pop[i+N1].k = 0;	
	}

	// set progress bar, initiate csv file, and start timer


	// simulation start
	struct virus next_gen[N0];
	struct virus *pop2;
	for (int repe=0; repe<rep; repe++){	
		for (int i=0; i<10; i++){ // run through generation
			mutate(pop);
			pop2 = step(pop,&next_gen);
			memcpy(pop,pop2,sizeof(struct virus)*N0); // cycle between pop and pop2 to continue looping.
			printf("pop0 outside of step id:%d, k1:%d, k2:%d, k:%d\n",pop[0].id,pop[0].k1,pop[0].k2,pop[0].k);
			printf("pop1 outside of step id:%d, k1:%d, k2:%d, k:%d\n",pop[1].id,pop[1].k1,pop[1].k2,pop[1].k);
			printf("pop2 outside of step id:%d, k1:%d, k2:%d, k:%d\n",pop[2].id,pop[2].k1,pop[2].k2,pop[2].k);
			printf("\n\n");
			//record
		}
	return 0;
}












void record(struct virus popop[]) {
	if (timestep) {

	} else {

	} 
}




void mutate(struct virus popop[]) {
	// goes through population and make every inidividual go through 
	// mutation process.
	// input: population struct array.
	// output: population struct with updated k values according to
	//   mutation rate.
	for (int i=0;i<N0; i++){
		int mut_num = bnldev(mu,L,&seed); // binomial pick of number of mutation based on mu.
		if (back == 1){ // back mutation
			if (popop[i].id == 1) { // individual is 1segment
				for (int i=0; i<mut_num; i++) { // add or subtract k
					if (ran1(&seed) < popop[i].k/L) { // k/L is the probability of back mutating.
						popop[i].k -= 1;
					} else {
						popop[i].k += 1;
					}
				}
			} else { // individual is 2segment
				for (int i=0; i<mut_num; i++) { // add or subtract k1 or k2, and k.
					if (ran1(&seed) < popop[i].k/L) {
						if (popop[i].k1==0) {
							popop[i].k2 -= 1;
							popop[i].k -= 1;
						} else if (popop[i].k2 == 0) {
							popop[i].k1 -= 1;
							popop[i].k -= 1;
						} else {
							if (ran1(&seed) < 0.5) {
								popop[i].k1 -= 1;
								popop[i].k -= 1;
							} else {
								popop[i].k2 -= 1;
								popop[i].k -= 1;
							}
						}
					} else {
						if (ran1(&seed) < 0.5) {
							popop[i].k1 += 1;
							popop[i].k += 1;
						} else {
							popop[i].k2 += 1;
							popop[i].k += 1;
						}
					}
				}
			}
		} else { // no back mutation 
			if (popop[i].id == 1) {
				popop[i].k += mut_num;
			} else {
				int breakpt = floor(ran1(&seed)*(mut_num+1));
				popop[i].k1 += breakpt;
				popop[i].k2 += mut_num - breakpt;
				popop[i].k = popop[i].k1 + popop[i].k2;
			}
		}
	}
}

struct virus *step(struct virus popop[],struct virus *next_gen_p) {
	// goes through reproduction process
	// the process is depicted at the top of the script.
	// input: pop struct array
	// output: next generation's pop struct array
	int l = 0; // next gen length
	while (l < N0) {
		int s1 = floor(ran1(&seed)*N0); // sample 1
		int s2 = floor(ran1(&seed)*N0); // sample 2
		//printf("s1:%d k1:%d k2:%d k:%d\n",s1,popop[s1].k1,popop[s1].k2,popop[s1].k);
		//printf("s2:%d k1:%d k2:%d k:%d\n",s2,popop[s2].k1,popop[s2].k2,popop[s2].k);
		if (popop[s1].id == 1 || popop[s2].id == 1) { // either parents is segment 1
			if (ran1(&seed) < 0.5) { // pick s1
				if (ran1(&seed) < pow(1.0-s,popop[s1].k)){
					next_gen_p[l].id = popop[s1].id;
					next_gen_p[l].k1 = popop[s1].k1;
					next_gen_p[l].k2 = popop[s1].k2;
					next_gen_p[l].k = popop[s1].k;
					l += 1;
				}
			} else { // pick s2
				if (ran1(&seed) < pow(1.0-s,popop[s2].k)){
					next_gen_p[l].id = popop[s2].id;
					next_gen_p[l].k1 = popop[s2].k1;
					next_gen_p[l].k2 = popop[s2].k2;
					next_gen_p[l].k = popop[s2].k;
					l += 1;
				}
			}
		} else { // both parents segment 2
			if (ran1(&seed) < 0.5) { // pick k1 from s1 and k2 from s2
				if (ran1(&seed) < pow(1.0-s,(popop[s1].k1 + popop[s2].k2))){
					next_gen_p[l].id = popop[s1].id;
					next_gen_p[l].k1 = popop[s1].k1;
					next_gen_p[l].k2 = popop[s2].k2;
					next_gen_p[l].k = popop[s1].k;
					l += 1;
				}
			} else { // pick k1 from s2 and k2 from s1
				if (ran1(&seed) < pow(1.0-s,(popop[s1].k2 + popop[s2].k1))){
					next_gen_p[l].id = popop[s2].id;
					next_gen_p[l].k1 = popop[s2].k1;
					next_gen_p[l].k2 = popop[s1].k2;
					next_gen_p[l].k = popop[s2].k;
					l += 1;
				}
			}
		}
	}
	printf("pop0 inside of step id:%d, k1:%d, k2:%d, k:%d\n",next_gen_p[0].id,next_gen_p[0].k1,next_gen_p[0].k2,next_gen_p[0].k);
	printf("pop1 inside of step id:%d, k1:%d, k2:%d, k:%d\n",next_gen_p[1].id,next_gen_p[1].k1,next_gen_p[1].k2,next_gen_p[1].k);
	printf("pop2 inside of step id:%d, k1:%d, k2:%d, k:%d\n",next_gen_p[2].id,next_gen_p[2].k1,next_gen_p[2].k2,next_gen_p[2].k);
	printf("\n\n");
	return next_gen_p;
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
