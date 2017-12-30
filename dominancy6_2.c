/* This model seeks for heterozygosity advantage in shifting seasons
through simulating population frequency each time step. 
depending on the expression level difference due to mutation allele, the 
survival probability changes along the normal curve.
In dominancy6, we have multiple mutations rising in a rate of
mutation in the loci. In addition the frequency now has stochastic
element through WF model simulation.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define _USE_MATH_DEFINES

float ran1(long *seed);
int intmax(int argc, int array[]);
int intmin(int argc,int array[]);
int intsum(int length, int array[]);
int *intseq(int init,int end, int inter);
double dnorm(double x, double mu, double sig);
double doublesum(int size, double a[]);

int genotype2index(int genotype[],int new_al_num){
  int index;
  if(intmin(2,genotype) == 1){
    index = intmax(2,genotype) - 1;
  } else {
    int *seq;
    seq = intseq(new_al_num,1,-1);
    int lesser = intmin(2,genotype);
    int *adds = (int *) malloc((lesser-1)*sizeof(int));
    for(int i=0; i<(lesser-1); i++){
      adds[i] = seq[i];
	}
	index = intsum(lesser-1,adds) + fabs(genotype[0]-genotype[1]);
  }
  return index;
}
int main(){
  FILE * fPointer;
  fPointer = fopen("result.csv","w");
  int N = 100; //population
  double mu = 1; //mutation rate
  double *x = (double *) malloc(sizeof(double));
  x[0] = 1.0; //frequency of alleles
  int *pop = (int *) malloc(sizeof(double));
  pop[0] = N; //array of population of genotypes
  //WRITE FILES FOR FREQ
  int t = 1; // generation amount
  int change = 5; //generation amount after a season change.
  double u1 = 0.333;
  double u2 = 0.666;
  double sig = 0.1665; //means of normal curves for different seasons.
  long seed = -1;
  double a1 = ran1(&seed)*0.5; //expression level of a wildtype
  double *all_exp = (double *) malloc(15*sizeof(double)); //expression level of each allele
  all_exp[0] = a1;
  int len_all_exp = 1;
  double *a = (double *) malloc(sizeof(double)); //expression level of all genotypes
  a[0] = 2*a1;
  int len_a = 1;
  int **genotypes_list = (int **) malloc(2*sizeof(int *));
  for(int i=0;i<2;i++){
  	genotypes_list[i] = (int *) malloc(sizeof(int));
  } 
  genotypes_list[0][0] = 1;
  genotypes_list[1][0] = 1;
  double *w1 = (double *) malloc(len_a*sizeof(double));
  double *w2 = (double *) malloc(len_a*sizeof(double));
  for(int i=0; i<len_a; i++){	
  	w1[i] = dnorm(a[i],u1,sig)/dnorm(u1,u1,sig);
  	w2[i] = dnorm(a[i],u2,sig)/dnorm(u2,u2,sig);
  }
  double *w = (double *) malloc(len_a*sizeof(double));
  memcpy(w,w1,len_a*sizeof(double));
  /*printf("a1 %f\n",a1);
  printf("all_exp %f\n",all_exp[0]);
  printf("a %f\n",a[0]);
  printf("genotypes_list %d%d\n",genotypes_list[0][0],genotypes_list[1][0]);
  printf("w1 %f\n",w1[0]);
  printf("w2 %f\n",w2[0]);
  printf("w %f\n",w[0]);*/
  //CHECK: IF MEMCPY FROM W2 OR W1 TO W WORKS PROPERLY
  for(int time=0; time<t; time++){
  	if(time%change == 1 && time>1){ //change season after specified generation time
  	  if(w[0] == w[1]){
  	  	free(w);
  	    memcpy(w,w2,len_a*sizeof(double));
  	  } else {
  	  	free(w);
  	    memcpy(w,w1,len_a*sizeof(double));
  	  }
  	}
    seed -= 1;
    float arise = ran1(&seed);
    //assumption: mutation arises only once per generation
    //protocol for when mutation arises
    if (arise < mu){ // algorithm for new mutation arising
      int *positive_pop_index = (int *) malloc(sizeof(int)*len_a); //select genotypes that have counts
      int j = 0;
      for(int i=0; i<len_a; i++){
        if(pop[i] > 0){
          positive_pop_index[j] = i;
          j += 1;
        }
	    }
      seed -= 1;
      int select_mutant = positive_pop_index[(int)floor(ran1(&seed)*(j+1))]; //select the genotype index that will mutate out of the genotypes with positive counts
      int mutant_arisen_genotype[2] = {genotypes_list[0][select_mutant],genotypes_list[1][select_mutant]}; //genotypes_list
      printf("mutant_arisen_genotype: %d %d\n",mutant_arisen_genotype[0],mutant_arisen_genotype[1]);
      seed -= 1;
      int one = (int)floor(ran1(&seed)*2);
      int arise_from = mutant_arisen_genotype[one]; //allele number that the mutation arose from
      seed -= 1;
      double d_new = ran1(&seed)-0.5; //new allele's expression dist from the original state
      double new_exp = all_exp[arise_from-1] + d_new; //new mutation expression level
      printf("d_new: %f, new_exp: %f\n", d_new, new_exp);
      if (new_exp<0){
      	new_exp = 0;
      }
      len_all_exp += 1; //length of allele expressions increased by 1
      all_exp[len_all_exp-1] = new_exp;
      printf("all_exp: %f %f %f\n",all_exp[0],all_exp[1],all_exp[2]);
      int num = intsum(len_all_exp,intseq(1,len_all_exp,1));
			for(int i=0; i<len_a; i++){
				free(genotypes_list[i]);
			}
			free(genotypes_list);
			free(a);
			a = (double *) malloc(num*sizeof(double));
			genotypes_list = (int **) malloc(2*sizeof(int *));
			for(int i=0; i<num; i++){
				genotypes_list[i] = (int *) malloc(num*sizeof(int));
			}
			int k=0;
			len_a = num; //update len_a
			for(int i=0; i<len_all_exp; i++){ //update genotypes' expression and genotype list
				for(int j=i; j<len_all_exp; j++){
					a[k] = all_exp[i] + all_exp[j];
					genotypes_list[0][k] = i+1;
					genotypes_list[1][k] = j+1;
					k++;
				}
			}
			printf("a:");
			for(int i=0; i<num; i++){
				printf("%f",a[i]);
			}
			printf("\n");
			printf("genotypes list:\n");
			for(int i=0; i<num; i++){
				for(int j=0; j<2; j++){
					printf("%d",genotypes_list[j][i]);
				}
				printf("\n");
			}
			if(w[0] == w[1]){
				free(w1);
				free(w2);
				free(w);
				w1 = (double *) malloc(len_a*sizeof(double));
				w2 = (double *) malloc(len_a*sizeof(double));
				w = (double *) malloc(len_a*sizeof(double));
				for(int i=0; i<len_a; i++){
					w1[i] = dnorm(a[i],u1,sig)/dnorm(u1,u1,sig);
					w2[i] = dnorm(a[i],u2,sig)/dnorm(u2,u2,sig);
				}
				memcpy(w,w1,len_a*sizeof(double));
			} else {
				free(w1);
				free(w2);
				free(w);
				w1 = (double *) malloc(len_a*sizeof(double));
				w2 = (double *) malloc(len_a*sizeof(double));
				w = (double *) malloc(len_a*sizeof(double));
				for(int i=0; i<len_a; i++){
					w1[i] = dnorm(a[i],u1,sig)/dnorm(u1,u1,sig);
					w2[i] = dnorm(a[i],u2,sig)/dnorm(u2,u2,sig);
				}
				memcpy(w,w1,len_a*sizeof(double));
			}
			x = (double *) realloc(x,len_all_exp*sizeof(double));
			x[arise_from-1] = x[arise_from-1] - (double)1/(2*N);
			x[len_all_exp-1] = (double)1/(2*N);
			printf("x: %f %f\n",x[0], x[1]);
			for(int i=0; i<len_a; i++){
				printf("pop[%d]= %d", i, pop[i]);
			}
			free(positive_pop_index);
    }
    double *x_square = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){ //get all the factors when x is squared
    	if(genotypes_list[0][i] == genotypes_list[1][i]){
    		x_square[i] = x[genotypes_list[0][i]]*x[genotypes_list[1][i]];
    	} else {
    		x_square[i] = 2*x[genotypes_list[0][i]]*x[genotypes_list[1][i]];
    	}
    }
    double *wx = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){
    	wx[i] = w[i]*x_square[i];
    }
    double wbar = doublesum(len_a,wx);
    free(pop);
    pop = (int *) malloc(len_a*sizeof(int));
    double *probs = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){
      probs[i] = wx[i]/wbar;
    }
    double *probs_accum = (double *) malloc(len_a*sizeof(double));
    probs_accum[0] = probs[0];
    for(int i=1; i<5; i++){
      probs_accum[i] = probs[i] + probs_accum[i-1];
    }
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
    free(probs);
    free(probs_accum);
    double factor_sum;
    for(int j=1; j<=len_all_exp; j++){
      factor_sum = 0;
      for(int i=0; i<len_a; i++){
        if(genotypes_list[0][i] ==j || genotypes_list[1][i] ==j){
          if(genotypes_list[1][i] == genotypes_list[0][i]){
            factor_sum += pop[i];
          } else {
            factor_sum += pop[i]/2;
          }
        }
      }
      x[j] = factor_sum/N;
    }
    for(int i=0; i<len_all_exp; i++){
      fprintf(fPonter,"%f,",x[i]);
    }
    fprintf("\n");
  }
  fclose(fPointer);
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

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++){
    sum += a[i];
  }
  return sum;
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

double dnorm(double x, double mu, double sig) {
  double y = 1/sqrt(2*M_PI*pow(sig,2))*exp(-1*pow((x-mu),2)/(2*pow(sig,2)));
  return y;
}

double doublesum(int size,double a[]){
  int i;
  float sum = 0;
  for(i=0; i<size; i++){
    sum += a[i];
  }
  return sum;
}