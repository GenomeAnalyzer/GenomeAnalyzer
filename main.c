#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.h"

int main(int argc, char const *argv[])
{
	//A : 00
	//T : 11
	//C : 10
	//G : 01
	//GTGAACTGGGCCCAT
	unsigned int sequence1[30] = {0,1,1,1,0,1,0,0,0,0,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0,0,0,1,1};

	//GTGAACTTGGACCAT
	unsigned int sequence2[30] = {0,1,1,1,0,1,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1};

	if(detecting_mutations(sequence1, 30))
		printf("Detected a probable mutation in sequence1\n");
	else
		printf("No mutation detected in sequence1\n");

	if(detecting_mutations(sequence2, 30))
		printf("Detected a probable mutation in sequence2\n");
	else
		printf("No mutation detected in sequence2\n");

	printf("ok\n");

	return 0;
}