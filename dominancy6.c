/* This model seeks for heterozygosity advantage in shifting seasons
through simulating population frequency each time step. 
depending on the expression level difference due to mutation allele, the 
survival probability changes along the normal curve.
In dominancy6, we have multiple mutations rising in a rate of
mutation in the loci. In addition the frequency now has stochastic
element through WF model simulation.*/

# include <stdio.h>

float ran1(long *seed);
int genotype2index(int geno[],int al_num);

int main(){
	int survivors[] = 
}