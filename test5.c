#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES

int main() {
	int **a = (int **) malloc(2*sizeof(int *));
	for(int i=0; i<2; i++){
		a[i] = (int *) malloc(2*sizeof(int));
	}

	for(int i=0; i<2; i++){
		free(a[i]); 
	}
	free(a);
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