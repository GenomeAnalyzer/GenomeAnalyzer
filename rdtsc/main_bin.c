#include <stdio.h>
#include <stdlib.h>
#include "rdtsc.h"
#include "../gene_bin.h"
#include <string.h>

#define MAX 31000
#define MAX_LINE 150
#define MAX_LOOP 500

void load_gene(char *filename, char *seq_char)
{
	char temp[MAX_LINE] = "";
	FILE *file;
	file = fopen(filename, "r");
	if(!file)
		printf("ERROR: load_gene: cannot open file %s\n", filename);

	if(file)
    {
		fgets(temp, MAX_LINE, file); 
        while(fgets(temp, MAX_LINE, file) != NULL)
        {
        	temp[strlen(temp)-1] = '\0';
        	strcat(seq_char, temp);
        }
    }
	//fscanf(file, "%s", seq_char);

	fclose(file);
}

int main(int argc, char *argv[])
{
	unsigned long long before = 0;
	unsigned long long after = 0;
	double elapsed = 0.0;

	char seq_char[MAX] = "";
	load_gene("LC528232.1.fasta", seq_char);
	// printf("%s", seq_char);
	unsigned long long seq_char_size = strlen(seq_char);

	char seq_char2[MAX] = "";
	load_gene("MN908947.3.fasta", seq_char2);
	// printf("%s", seq_char2);
	unsigned long long seq_char_size2 = strlen(seq_char2);


	// Variable used for function
	long int *seq_long = NULL;
	long int *seq_long2 = NULL;
	float cms = 0;
	char *rna_seq_long = NULL;
	char *aa_seq_long = NULL;
	
	seq_long2 = convert_to_binary(seq_char2, seq_char_size2);
  	
  	gene_map_t g;
    g.gene_start = malloc(sizeof(*g.gene_start) * seq_char_size * 2);
    g.gene_end = malloc(sizeof(*g.gene_end) * seq_char_size * 2);
    
    mutation_map m;
	m.size = malloc(sizeof(unsigned long) * 5);
	m.start_mut = malloc(sizeof(unsigned long) * 5);
	m.end_mut = malloc(sizeof(unsigned long) * 5);
	for(int i = 0; i < 5; i++)
	{
		m.size[i]=0;
    	m.start_mut[i]=0;
    	m.end_mut[i]=0;   	
    }

    printf("Binaries Functions\t    | Cycles\n");
    printf("-----------------------------------------\n");

	/*-----convert_to_binary-----*/
	for(int i = 0; i < MAX_LOOP; i++)
	{
    	before = rdtsc();
    	seq_long = set_binary_array(seq_char, seq_char_size);
    	after = rdtsc();

    	elapsed += (double)(after - before);
	}
	printf("convert_to_binary\t    : %.3lf\n", elapsed / MAX_LOOP);
	elapsed = 0;

	// /*-----generating_mRNA-----*/
	for(int i = 0; i < MAX_LOOP; i++)
	{
    	before = rdtsc();
		rna_seq_long = generating_mRNA(seq_long, 0, seq_char_size);
    	after = rdtsc();

    	elapsed += (double)(after - before);
	}
	printf("generating_mRNA\t\t    : %.3lf\n", elapsed / MAX_LOOP);
	elapsed = 0;
	
	// /*-----detecting_genes-----*/
	// for(int i = 0; i < MAX_LOOP; i++)
	// {
	// 	before = rdtsc();
	// 	detecting_genes(seq_long, seq_char_size, &g);
	// 	after = rdtsc();

	// 	elapsed += (double)(after - before);
	// }
	// printf("detecting_genes : %lf cycles\n", elapsed / MAX_LOOP);
	// elapsed = 0;
	
	// /*-----generating_amino_acid_chain-----*/
	for(int i = 0; i < MAX_LOOP; i++)
	{
    	before = rdtsc();
		aa_seq_long = generating_amino_acid_chain(seq_long, 0, seq_char_size);
    	after = rdtsc();

    	elapsed += (double)(after - before);
	}
	printf("generating_amino_acid_chain : %.3lf\n", elapsed / MAX_LOOP);
	elapsed = 0;

	// /*-----detecting_mutations-----*/
	for(int i = 0; i < MAX_LOOP; i++)
	{
    	before = rdtsc();
    	detecting_mutations(seq_long, 0, seq_char_size, m);
    	after = rdtsc();

    	elapsed += (double)(after - before);
	}
	printf("detecting_mutations\t    : %.3lf\n", elapsed / MAX_LOOP);
	elapsed = 0;

	// /*-----calculating_matching_score-----*/
	for(int i = 0; i < MAX_LOOP; i++)
	{
    	before = rdtsc();
		cms = calculating_matching_score(seq_long, 0, seq_char_size, seq_long2, 0, 2 * seq_char_size2);
    	after = rdtsc();

    	elapsed += (double)(after - before);
	}
	printf("calculating_matching_score  : %.3lf\n", elapsed / MAX_LOOP);
	elapsed = 0;

	printf("\n");

	// free
	free(g.gene_start);
	free(g.gene_end);
	free(m.size);
	free(m.start_mut);
	free(m.end_mut);

	return 0;
}
