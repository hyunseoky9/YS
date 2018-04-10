#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int back = 0;
int timestep = 0;
int krecord = 0;
int rep = 1;
int L = 300;
double s = 0.05;
#define N0 1000
int K = 1000;
float mu = 0.0005;
int gen_num = 10;
double cost = 0.00;
double r = 0.5;
double N1r = 0.5;
long seed = 0;
long bnlseed = 0;

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
	struct virus pop[N0];
	//initialize
	int N1 = N0*N1r;
	int N2 = N0*(1-N1r);
	for (int i=0;i<N1;i++){
		pop[i].id = 1;
		pop[i].k1 = 0;
		pop[i].k2 = 0;
		pop[i].k = 0;
	}
	for (int i=0;i<N2;i++){
		pop[i+500].id = 2;
		pop[i+500].k1 = 0;
		pop[i+500].k2 = 0;
		pop[i+500].k = 0;	
	}
	printf("pop1 has %d segments\n",pop[0].id);
	printf("pop500 has %d segments\n",pop[500].id);

	//go through generations
	struct virus next_gen[N0];
	struct virus *pop2;

	//mutate
	pop2 = step(pop,&next_gen);
	pop2[N0-1].id = 500;
	printf("pop2 segment: %d",pop2[N0-1].id);
	memcpy(pop,pop2,sizeof(struct virus)*N0);

	pop2 = step(pop,&next_gen);

	return 0;
}


void mutate(struct virus popop[]) {
	for (int i=0;i<N0; i++){
		bnlseed -= 1;
		int mut_num = bnldev(mu,L,&bnlseed);
		if (back == 1){ // back mutation

		} else { // no back mutation 
			if (popop[i].id == 1) {
				popop[i].k += mut_num;
			} else {
				seed += 1;
				int breakpt = floor(ran1(&seed)*(mut_num+1));
				popop[i].k1 += breakpt;
				popop[i].k2 += mut_num - breakpt;
				popop[i].k = 
			}
		}
	}
}


struct virus *step(struct virus popop[],struct virus *next_gen_p) {
	int l = 0; // next gen length
	for (int i=0;i<N0;i++) {
		next_gen_p[i].id = popop[i].id;
		next_gen_p[i].k1 = popop[i].k1;
		next_gen_p[i].k2 = popop[i].k2;
		next_gen_p[i].k = popop[i].k;
	}
	/*while (l <  N0) {
		seed += 1;
		s1 = (int)floor(ran1(&seed)); // sample 1
		seed += 1;
		s2 = (int)floor(ran1(&seed)); // sample 2
		if (pop[s1].id == 1 || pop[s2].id == 1) {
			seed += 1;
			if (ran1(&seed) < 0.5){ // pick sample 0 or 1
				seed += 1;
				if (ran1(&seed) < pow(1.0-s,pop[s1].k)){
						next_gen[l].id = pop[s1].id;
						next_gen[l].k1 = pop[s1].k1;
						next_gen[l].k2 = pop[s1].k2;
						next_gen[l].k = pop[s1].k;
				}
			}
		}
	}*/
	printf("pointer in function %p\n", next_gen_p);
	return next_gen_p;
}



/*

		if sample[0].id == 1 or sample[1].id == 1:
			if np.random.uniform(0,1) < 0.5: # sample 0 or 1
				if np.random.uniform(0,1) < (1-s)**sample[0].k: # progeny live or die
					if sample[0].id == 1:
						next_gen.append(Virus1(sample[0].k))
						popsize[0] += 1
					else:
						next_gen.append(Virus2(sample[0].k1,sample[0].k2))
						popsize[1] += 1
			else:
				if np.random.uniform(0,1) < (1-s)**sample[1].k:
					if sample[1].id == 1:
						next_gen.append(Virus1(sample[1].k))
						popsize[0] += 1
					else:
						next_gen.append(Virus2(sample[1].k1,sample[1].k2))
						popsize[1] += 1
		elif sample[0].id == 2 and sample[1].id == 2:
			if np.random.uniform(0,1) < 0.5:
				if np.random.uniform(0,1) < (1-s)**(sample[0].k1+sample[1].k2):
					next_gen.append(Virus2(sample[0].k1,sample[1].k2))
					popsize[1] += 1
			else:
				if np.random.uniform(0,1) < (1-s)**(sample[0].k2+sample[1].k1):
					next_gen.append(Virus2(sample[0].k2,sample[1].k1))
					popsize[1] += 1
	print('sampling #:',samplenum)
	if np.sum(popsize) != N0:
		raise ValueError('popsize doesn\'t add up to N0!! sum=%s'%(next_gen))
	return next_gen, popsize


*/


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
