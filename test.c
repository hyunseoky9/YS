#include <stdio.h>
#include <stdlib.h>
int main(void)
{
    // This program will create same sequence of 
    // random numbers on every program run 
 
    for(int i = 0; i<5; i++)
        printf(" %f ", (float)rand()/RAND_MAX);

    return 0;
}	 	
