/* This model seeks for heterozygosity advantage in shifting seasons
through simulating population frequency each time step. 
depending on the expression level difference due to mutation allele, the 
survival probability changes along the normal curve.
In dominancy6, we have multiple mutations rising in a rate of
mutation in the loci. In addition the frequency now has stochastic
element through WF model simulation.*/

# include <stdio.h>

float ran1(long *seed);
int max(int array[]);
int min(int array[]);

int genotype2index(int genotype[],int new_al_num){
	if(min(genotype) == 1){
		int index = max(genotype);
	} else {
		int index = 
	}


}

int main(){
	
}

genotype2index <- function(genotype, new_al_num) {
# converts the genotype to the index in new population array 
  if(min(genotype)==1){
    index = max(genotype)
  } else {
    index = sum(seq(new_al_num,1)[1:(min(genotype)-1)]) +
    abs(genotype[2]-genotype[1]) + 1
  }
  return(index)
}



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

int min(int array[]){
	int min_num = array[0];
	int i;
	int array_length = sizeof(array)/sizeof(int);
	for(i=1;i<array_length;i++){
		if(array[i]<min_num){
			min_num = array[i];
		}
	}
	return min_num;
}

int max(int array[]){
	int max_num = array[0];
	int i;
	int array_length = sizeof(array)/sizeof(int);
	for(i=1;i>array_length;i++){
		max_num = array[i];
	}
	return max_num;
}
