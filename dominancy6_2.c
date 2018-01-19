/* This model seeks for polymorphisms caused by shifting seasons.
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
#define PRINTOUT 0

float ran1(long *seed); //random number generator between 0 and 1.
int intmax(int argc, int array[]); //max value of an integer array
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int length, int array[]); //sum of an integer array
int *intseq(int init,int end, int inter); //sequence of integer given by initial, end, and interval values.
double dnorm(double x, double mu, double sig); //value of a norml pdf function given a parameter.
double doublesum(int size, double a[]); //sum of a double array

int main(){

  //define the parameters.
  FILE * fPointer;
  fPointer = fopen("result.csv","w"); //file output of frequency in each time step.
  int N = 1000; //population
  double mu = 0.00001; //mutation rate per generation
  double *x = (double *) malloc(sizeof(double)); //allele frequency
  x[0] = 1.0; //frequency of alleles
  int *pop = (int *) malloc(sizeof(int));
  pop[0] = N; //array of population of genotypes
  int t = 10000; // generation amount
  int change = 5; //generation amount after a season change.
  double u1 = 0.333; //mean of the normal curve for season1.
  double u2 = 0.666; //mean of the normal curve for season2.
  double sig = 0.1665; //sd for normal curves for both seasons.
  long seed;
  printf("ENTER SEED VALUE (negative integer): ");
  scanf("%ld", &seed);
  printf("\n %ld",seed);
  double a1 = ran1(&seed)*0.5; //expression level of a wildtype allele
  //double a1 = 0.5; //used when there's NO SELECTION
  double *all_exp = (double *) malloc(15*sizeof(double)); //expression level of each allele
  all_exp[0] = a1;
  int len_all_exp = 1;
  double *a = (double *) malloc(sizeof(double)); //expression level of all genotypes
  a[0] = 2*a1;
  int len_a = 1; //length of a
  int **genotypes_list = (int **) malloc(2*sizeof(int *)); //list of all the genotypes at current time step each column is a genotype.
  for(int i=0;i<2;i++){
  	genotypes_list[i] = (int *) malloc(sizeof(int));
  }
  genotypes_list[0][0] = 1;
  genotypes_list[1][0] = 1;
  double *w1 = (double *) malloc(len_a*sizeof(double)); //fitness of a genotype at season1
  double *w2 = (double *) malloc(len_a*sizeof(double)); //fiteness of a genotype at season2
  for(int i=0; i<len_a; i++){	
  	w1[i] = dnorm(a[i],u1,sig)/dnorm(u1,u1,sig);
  	w2[i] = dnorm(a[i],u2,sig)/dnorm(u2,u2,sig);
  }
  double *w = (double *) malloc(len_a*sizeof(double)); //current fitness of a genotype
  memcpy(w,w1,len_a*sizeof(double)); //current fitness starts with season1's fitness

  //Loop each time step
  for(int time=0; time<t; time++){
    if(PRINTOUT){printf("\n\n\n -----------------------time=%d------------------------\n\n",time);}
  	
    //change season after specified generation time 'change'
    if((time+1)%change == 1 && (time+1)>1){
      if(PRINTOUT){printf("TIME TO CHANGEEEEE\n");}
  	  if(w[0] == w1[0]){
        if(PRINTOUT){printf("CHANGE BACK TO SEASON2\n");}
  	  	//free(w);
  	    memcpy(w,w2,len_a*sizeof(double));
  	  } else {
        if(PRINTOUT){printf("CHANGE BACK TO SEASON1\n");}
  	  	//free(w);
  	    memcpy(w,w1,len_a*sizeof(double));
  	  }
  	}

    //steps when mutation arises after a time step. 
    //allele that is to mutate is chosen, and frequency 'x', genotype_list, a, and all_exp is updated.
    //assumption: mutation arises only once in a time step if it occurs
    seed -= 1;
    float arise = ran1(&seed);
    if (arise < 2*N*mu){ 
      
      //select the gene to mutate
      int *positive_pop_index = (int *) malloc(len_a*sizeof(int)); //select genotypes that have positive counts
      int j = 0;
      for(int i=0; i<len_a; i++){
        if(pop[i] > 0){
          positive_pop_index[j] = i;
          j += 1;
        }
	    }
      seed -= 1;
      int select_mutant = positive_pop_index[(int)floor(ran1(&seed)*j)]; //select the genotype index that will mutate out of the genotypes with positive counts
      int mutant_arisen_genotype[2] = {genotypes_list[0][select_mutant],genotypes_list[1][select_mutant]}; //genotypes_list
      if(PRINTOUT){printf("mutant_arisen_genotype: %d %d\n",mutant_arisen_genotype[0],mutant_arisen_genotype[1]);}
      seed -= 1;
      int one = (int)floor(ran1(&seed)*2);
      int arise_from = mutant_arisen_genotype[one]; //allele index that the mutation arose from
      
      //give an expression level for the new allele
      seed -= 1;
      double d_new = ran1(&seed)-0.5; //new allele's expression dist from the original allele's
      double new_exp = all_exp[arise_from-1] + d_new; //new mutation expression level
      //double new_exp = 0.5; //allele expression level when there's NO SELECTION
      if(PRINTOUT){printf("d_new: %f, new_exp: %f\n", d_new, new_exp);}
      if (new_exp<0){
      	new_exp = 0; 
      }
      len_all_exp += 1; 
      all_exp[len_all_exp-1] = new_exp; 

      //update genotype_list and genotype expression array.
      int *sequence = intseq(1,len_all_exp,1); 
      int num = intsum(len_all_exp, sequence); 
      free(sequence); 
			for(int i=0; i<2; i++){
				free(genotypes_list[i]);
			}
			free(genotypes_list);
			free(a);
			a = (double *) malloc(num*sizeof(double));
			genotypes_list = (int **) malloc(2*sizeof(int *));
			for(int i=0; i<2; i++){
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
      if(PRINTOUT){
  			printf("a:");
  			for(int i=0; i<len_a; i++){
  				printf("%f ",a[i]);
  			}
  			printf("\n");
  			printf("genotypes list:\n");
  			for(int i=0; i<len_a; i++){
  				for(int j=0; j<2; j++){
  					printf("%d",genotypes_list[j][i]);
  				}
  				printf("\n");
  			}
      }
      //update genotype fitness array
      if(w[0] == w1[0]){
        if(PRINTOUT){printf("W IS W1\n");}
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
        if(PRINTOUT){printf("W IS W2\n");}
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
				memcpy(w,w2,len_a*sizeof(double));
			}

      //update freqeuncy array
			x = (double *) realloc(x, len_all_exp*sizeof(double));
			x[arise_from-1] = x[arise_from-1] - (double)1/(2*N);
			x[len_all_exp-1] = (double)1/(2*N);
      if(PRINTOUT){
        printf("x: ");
        for(int i=0; i<len_all_exp;i++){
          printf("%f ",x[i]);
        }
        printf("\n");
      }
			free(positive_pop_index);
    }

    //calculate genotype ratio
    double *x_square = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){ //get all the factors when x is squared
    	if(genotypes_list[0][i] == genotypes_list[1][i]){
    		x_square[i] = x[genotypes_list[0][i]-1]*x[genotypes_list[1][i]-1];
    	} else {
    		x_square[i] = 2*x[genotypes_list[0][i]-1]*x[genotypes_list[1][i]-1];
    	}
    }
    double *wx = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){
    	wx[i] = w[i]*x_square[i];
    }
    if(PRINTOUT){
      printf("x_square: ");
      for(int i=0; i<len_a; i++){
        printf("%f ",x_square[i]);
      }
      printf("\n");
      printf("w: ");
      for(int i=0; i<len_a; i++){
        printf("%f ",w[i]);
      }
      printf("\n");
    }
    free(x_square);
    double wbar = doublesum(len_a,wx);

    //fill in the population array for the next generation based on the probability from the genotype ratio
    free(pop);
    pop = (int *) malloc(len_a*sizeof(int));
    for(int i=0; i<len_a; i++){
      pop[i] = 0;
    }
    double *probs = (double *) malloc(len_a*sizeof(double));
    for(int i=0; i<len_a; i++){
      probs[i] = wx[i]/wbar;
    }
    if(PRINTOUT){
      printf("probs: ");
      for(int i=0; i<len_a; i++){
        printf("%f ",probs[i]);
      }
      printf("\n");
    }
    free(wx);
    double *probs_accum = (double *) malloc(len_a*sizeof(double));
    probs_accum[0] = probs[0];
    for(int i=1; i<len_a; i++){
      probs_accum[i] = probs[i] + probs_accum[i-1];
    }
    if(PRINTOUT){
      printf("probs_accum: ");
      for(int i=0; i<len_a; i++){
        printf("%f ",probs_accum[i]);
      }
      printf("\n");
    }
    for(int i=0; i<N; i++){
      seed -= 1;
      float val = ran1(&seed);
      if(val<probs_accum[0]){
        pop[0] += 1;
      } else {
        for(int j=1; j<len_a; j++){
          if(val<probs_accum[j] && val>probs_accum[j-1]){
            pop[j] += 1;
          }
        }
      }
    }
    if(PRINTOUT){
      printf("pop: ");
      for(int i=0; i<len_a; i++){
        printf("%d ",pop[i]);
      }
      printf("\n");
    }
    free(probs);
    free(probs_accum);

    //update frequency array based on the new population array
    double factor_sum;
    for(int j=1; j<=len_all_exp; j++){
      factor_sum = 0;
      for(int i=0; i<len_a; i++){
        if(genotypes_list[0][i] == j || genotypes_list[1][i] == j){
          if(genotypes_list[0][i] == genotypes_list[1][i]){
            factor_sum += pop[i];
          } else {
            factor_sum += pop[i]/(double)2;
          }
        }
      }
      x[j-1] = factor_sum/N;
    }
    if(PRINTOUT){
      printf("x: ");
      for(int i=0; i<len_all_exp; i++){
        printf("%f ",x[i]);
      }
      printf("\n");
    }
    for(int i=0; i<len_all_exp; i++){ //write out the frequency array
      fprintf(fPointer,"%.3f,",x[i]); 
    }
    fprintf(fPointer,"\n");
  }
  free(x);
  free(pop);
  free(all_exp);
  free(a);
  for(int i=0; i<2; i++){
    free(genotypes_list[i]);
  }
  free(genotypes_list);
  free(w1);
  free(w2);
  free(w);
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
#undef  NTAB
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