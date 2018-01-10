#include <stdio.h>
#include <stdlib.h>

int main() {
	FILE * fPointer;
	fPointer = fopen("testtesttest.csv","w");
	fprintf(fPointer,"test %.2f",0.123f);
	fclose(fPointer);
	return 0;
}