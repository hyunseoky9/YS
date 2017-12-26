#include <stdio.h>
#include <math.h>
int intseq(int init, int end, int inter){
	double a = fabs(end-init)/inter;
	int add_amount = floor(a);
	int sequence[add_amount+1];
	sequence[0] = init;
	if(add_amount>0){x
		int i;
		for(i=1;i<add_amount;i++){
			sequence[i] = sequence[i] + inter;
		}
	}
	return sequence;
}

int main(){
	int *seq = intseq(2,6,2);
	printf("sequence goes %d %d %d",seq[0],seq[1],seq[2]);
	return 0;
}