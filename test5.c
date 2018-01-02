#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES

int main() {
	for(int i=0; i<4; i++){
		double *x = (double *) malloc(sizeof(double));
		x[0] = 0.1;
		x = (double *) realloc(x,2*sizeof(double));
		free(x);
		printf("yep");
	}
	/*a = (int **) malloc(2*sizeof(int *));
	for(int i=0; i<2; i++){
		a[i] = (int *) malloc(2*sizeof(int));
	}*/

	/*for(int i=0; i<2; i++){
		free(a[i]); 
	}
	free(a);*/

	return 0;
}