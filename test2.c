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
	int len_all_exp = 4;
	int len_a = 6;
	int num = intsum(len_all_exp,intseq(1,len_all_exp,1));
	int *newpop = (int *) malloc(num*sizeof(int));
	int j = 0;
	int pop[6] = {1,2,3,4,5,6};
	int *tempop = (int *) malloc((len_a)*sizeof(int));
	int len_tempop = len_a;
	memcpy(tempop,pop,len_a*sizeof(int));
	for(int i=(len_all_exp-1); i>0; i--){
      	for(int k=0;k<i;k++){
      		newpop[j+k] = tempop[k];
        }
        int *temp = (int *) malloc((len_tempop-i)*sizeof(int));
        //printf("temp: %d\n",temp[0]);
        int l = 0;
        for(int k=i;k<len_tempop;k++){
        	temp[l] = tempop[k];
        	l += 1;
        	//printf("temp: %d\n",temp[0]);
        	//printf("k is %d\n",k);
        }
        len_tempop = l;
        free(tempop);
        tempop = (int *) malloc(len_tempop*sizeof(int));        
        if(l>0){
        	memcpy(tempop, temp, len_tempop*sizeof(int));
    	}
    	free(temp);
    	//printf("tempop: %d %d %d %d %d %d\n",tempop[0],tempop[1],tempop[2],tempop[3],tempop[4],tempop[5]);
	  	j = j + (i-1) +2; 
	}
	printf("newpop\n");
	for(int i=0; i<num; i++){
		printf("%d\n",newpop[i]);
	}
	FILE * fPointer;
	fPointer = fopen("test2.csv","w");
	fprintf(fPointer,"%f,%f\n",3.4f,3.44f);
	fprintf(fPointer,"%f,%f\n",0.654f,0.235f);
	fclose(fPointer);
	//printf("index is %d", index);
	return 0;
}