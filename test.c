#include <stdio.h>
#include <stdlib.h>

int main(){
	FILE *ofp;
	ofp = fopen("fuck it.txt","w");
	fprintf(ofp,"Just joking, we don't\n");
	fprintf(ofp,"new line\n");
	fclose(ofp);
	return 0;
}