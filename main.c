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
	unsigned short sequence1[30] = {0,1,1,1,0,1,0,0,0,0,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0,0,0,1,1};

	//GTGAACTTGGACCAT
	unsigned short sequence2[30] = {0,1,1,1,0,1,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1};

	mutation_map M;
	M.size = malloc(30 * sizeof(unsigned long));
	M.start_mut = malloc(30 * sizeof(unsigned long));
	M.end_mut = malloc(30 * sizeof(unsigned long));

	detecting_mutations(sequence1, 30, M);
	if(M.size[0]){
		printf("Detected a probable mutation in sequence1\n");
		printf("size: %ld, start: %ld, end: %ld\n", M.size[0], M.start_mut[0], M.end_mut[0]);
	}
	else
		printf("No mutation detected in sequence1\n");

	M.size[0] = 0;
    M.end_mut[0] = 0;
    M.size[0] = 0;
	detecting_mutations(sequence2, 30, M);
	if(M.size[0]){
		printf("Detected a probable mutation in sequence2\n");
		printf("size: %ld, start: %ld, end: %ld\n", M.size[0], M.start_mut[0], M.end_mut[0]);
	}
	else
		printf("No mutation detected in sequence2\n");

	printf("ok\n");

	return 0;
}