#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES

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
	int N = 100;
	double probs[5] = {0.1,0.3,0.2,0.2,0.2};
	double *probs_accum = (double *) malloc(5*sizeof(double));
	probs_accum[0] = probs[0];
	for(int i=0; i<5; i++){
		probs_accum[i] += probs_accum[i-1];
	}
	for(int i=0,i<N; )
	//printf("index is %d", index);
	return 0;
}