#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES

int main() {
	int *a = (int *) malloc(sizeof(int));
	int *b = (int *) malloc(sizeof(int));
	a[0] = 2;
	b[0] = 2;
	//free(a);
	//free(a);
	free(a);
	free(b);
	printf("%d", b[0]);
	return 0;
}