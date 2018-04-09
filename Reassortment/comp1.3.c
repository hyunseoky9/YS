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
int mu = 0.0005;
int gen_num = 10;
double cost = 0.00;
double r = 0.5;
double N1r = 0.5;
long seed = 0;

struct virus {
	int id;
	int k1;
	int k2;
	int k;
};

//https://stackoverflow.com/questions/2620146/how-do-i-return-multiple-values-from-a-function-in-c
float ran1(long *seed);
struct virus *step(struct virus popop[],struct virus *next_gen_p);

int main(void) {
	struct virus pop[100];
	pop[0].id = 1;
	pop[0].k1 = 2;
	pop[0].k2 = 0;
	pop[0].k = 2;
	struct virus next_gen[N0];
	printf("pop0 is %d segmented and have %d k1 and %d k mutations\n",pop[0].id,pop[0].k1,pop[0].k);
	struct virus *pop2;
	pop2 = step(pop,&next_gen);
	printf("pointer after function %p\n", pop2);
	printf("Now pop0 is %d segmented and have %d k1 and %d k mutations\n",pop2[0].id,pop2[0].k1,pop2[0].k);
	return 0;
}

struct virus *step(struct virus popop[],struct virus *next_gen_p) {
	int l = 0; // declared next gen length
	next_gen_p[0].id = popop[0].id;
	next_gen_p[0].k1 = popop[0].k1;
	next_gen_p[0].k2 = popop[0].k2;
	next_gen_p[0].k = popop[0].k;
	/*while (len_next_gen <  N0) {
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
