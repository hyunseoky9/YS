#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES

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

int intmax(int argc, int array[]){
  int max_num = array[0];
  int i;
  for(i=1;i<argc;i++){
  	if(array[i] > max_num){
    	max_num = array[i];
    }
  }
  return max_num;
}

int *intseq(int init, int end, int inter){
	double a = (double)fabs(end-init)/fabs(inter);
	int add_amount = floor(a);
	int *sequence = (int *) malloc((add_amount+1)*sizeof(int));
	sequence[0] = init;
	if(add_amount>0){
		int i;
		for(i=1;i<(add_amount+1);i++){
			sequence[i] = sequence[i-1] + inter;
		}
	}
	return sequence;
}

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++){
    sum += a[i];
  }
  return sum;
}

int genotype2index(int genotype[],int new_al_num){
	int index;
	if(intmin(2,genotype) == 1){
		index = intmax(2,genotype);
	} else {
		int *seq;
		seq = intseq(new_al_num,1,-1);
		int lesser = intmin(2,genotype);
		int *adds = (int *) malloc((lesser-1)*sizeof(int));
		for(int i=0; i<(lesser-1); i++){
			adds[i] = seq[i];
		}
		index = intsum(lesser-1,adds) + fabs(genotype[0]-genotype[1]) + 1;
	}
	return index;
}

double dnorm(double x, double mu, double sig) {
	double y = 1/sqrt(2*M_PI*pow(sig,2))*exp(-1*pow((x-mu),2)/(2*pow(sig,2)));
	return y;
}

int main(){
	long seed = -1;
	int N = 100;
	double probs[5] = {0.1,0.3,0.2,0.2,0.2};
	int pop[5] = {0,0,0,0,0};
	double *probs_accum = (double *) malloc(5*sizeof(double));
	probs_accum[0] = probs[0];
	for(int i=1; i<5; i++){
		probs_accum[i] = probs[i] + probs_accum[i-1];
	}
	printf("probs_accum:");
	for(int  i=0; i<5; i++){
		printf("%f ", probs_accum[i]);
	}
	printf("\n");
	for(int i=0; i<N; i++){
		seed -= 1;
		if(ran1(&seed)<probs_accum[0]){
			pop[0] += 1;
		} else {
			for(int j=1; j<5; j++){
				if(ran1(&seed)<probs_accum[j] && ran1(&seed)>probs_accum[j-1]){
					pop[j] += 1;
				}
			}	
		}
	}
	printf("Pop:");
	for(int i=0; i<5; i++){
		printf("%d ", pop[i]);
	}
	return 0;
}