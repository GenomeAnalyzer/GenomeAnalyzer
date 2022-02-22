#include <nmmintrin.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "versions/bin/headers/gene_bin.h"
#include "rdtsc/rdtsc.h"

#define MAX 36
#define MAX_LINE 150
#define MAX_ADN 100000

/*
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
    fclose(file);
}
*/

unsigned long long load_from_file(char *filename, char *seq_char)
{
	FILE *file;
	file = fopen(filename, "r");
	if(!file)
		printf("ERROR: load_from_file: cannot open file 'seq.txt'\n");

	unsigned long long i = 0;
	while((seq_char[i] = fgetc(file)) != EOF)
		i++;
	printf("size seq_char  = %llu\n", i-1);

	fclose(file);
	return i-1;
}

#pragma GCC target("avx2")
void detecting_genes_2(const long int *gene, const long int gene_size, gene_map_t* gene_map) {
        gene_map->genes_counter = 0;

    // Check if memory ever have been allocated and allocate it if not
    if(!gene_map->gene_start || !gene_map->gene_end){
        gene_map->gene_start = malloc(sizeof(*gene_map->gene_start) * MAX_GENES);
        gene_map->gene_end = malloc(sizeof(*gene_map->gene_end) * MAX_GENES);

        if (!gene_map->gene_start || !gene_map->gene_end){
            printf("ERROR: detecting_genes: cannot allocate memory\n");
            return;
        }
    }

    // Start codon
    __m128i AUG = _mm_set_epi16(0, 0, 1, 1, 0, 1, 0, 0);
    // Stop codons
    __m128i UAA = _mm_set_epi16(1, 1, 0, 0, 0, 0, 0, 0);
    __m128i UGA = _mm_set_epi16(1, 1, 0, 1, 0, 0, 0, 0);
    __m128i UAG = _mm_set_epi16(1, 1, 0, 0, 0, 1, 0, 0);

    __m128i codon_to_test = _mm_setzero_si128();

    int start_pos = -1;
    long int i = 0;
   //Parse the binary array, and find all the start and stop codons
    while ((i + 6) <= gene_size) {
        // Each nucleotides can be A, U, G or C
        codon_to_test = _mm_set_epi16(get_binary_value(gene, i),
                                        get_binary_value(gene, i + 1),
                                        get_binary_value(gene, i + 2),
                                        get_binary_value(gene, i + 3),
                                        get_binary_value(gene, i + 4),
                                        get_binary_value(gene, i + 5),
                                        0, 0);

        __m128i cmp_AUG = _mm_xor_si128(codon_to_test, AUG);

	__m128i cmp_UAA = _mm_xor_si128(codon_to_test, UAA);
        __m128i cmp_UAG = _mm_xor_si128(codon_to_test, UAG);
        __m128i cmp_UGA = _mm_xor_si128(codon_to_test, UGA);

        //If a start pos and a stop pos doesn't exist, search for AUG
        if (_mm_testz_si128(cmp_AUG, cmp_AUG)) {
        //if AUG, it's the start of a gene
            start_pos = i;
            i += 6;
        }
        else if ((start_pos != -1 ) 
	         && (_mm_testz_si128(cmp_UAA, cmp_UAA)
	             || _mm_testz_si128(cmp_UAG, cmp_UAG)
                     || _mm_testz_si128(cmp_UGA, cmp_UGA))) {
             //It's the end of a gene          
             //If a start pos and an stop pos has been found, a gene exists so we save it in the struc
            gene_map->gene_start[gene_map->genes_counter] = start_pos;
            gene_map->gene_end[gene_map->genes_counter] = i+5;

	    gene_map->genes_counter++;

	    start_pos = -1;
	    i += 6;
	}
        else
            i += 2;
    }
}



#pragma GCC target("avx2")
void test_m128i()
{
    int test = 0;
    __m128i v0 = _mm_setzero_si128();
    __m128i v1 = _mm_setzero_si128();
    __m128i vcmp = _mm_setzero_si128();

    v0 = _mm_set_epi16(0,1,0,1,0,1,1,1);
    v1 = _mm_set_epi16(1,1,1,1,1,1,0,0);

    vcmp = _mm_xor_si128(v0, v1);
    test = _mm_testz_si128(vcmp, vcmp);

    if(test)
        printf("equal: %d\n\n", test);
    else
        printf("different: %d\n\n", test);
}

unsigned long long test1(unsigned long long seq_char_size, long int *seq_long)
{
    gene_map_t g;
    g.gene_start = malloc(sizeof(*g.gene_start) * seq_char_size * int_SIZE);
    g.gene_end = malloc(sizeof(*g.gene_end) * seq_char_size * int_SIZE);

    detecting_genes(seq_long, 2 * seq_char_size, &g);
    //printf("%llu genes\n\n", g.genes_counter);

    free(g.gene_start);
    free(g.gene_end);
    return g.genes_counter;
}

unsigned long long test2(unsigned long long seq_char_size, long int *seq_long)
{
    gene_map_t g;
    g.gene_start = malloc(sizeof(*g.gene_start) * seq_char_size * int_SIZE);
    g.gene_end = malloc(sizeof(*g.gene_end) * seq_char_size * int_SIZE);

    detecting_genes_2(seq_long, 2 * seq_char_size, &g);
    //printf("%llu genes\n\n", g.genes_counter);

    free(g.gene_start);
    free(g.gene_end);
    return g.genes_counter;
}

int main(int argc, char *argv[])
{
    // --- TEST AVX2 --- //
    test_m128i();

    // --- TESTS DETECTING GENES --- //
    char *seq_char = NULL;
    unsigned long long seq_char_size = 0;
    long int *seq_long = NULL;

    seq_char = malloc(sizeof(*seq_char) * MAX_ADN);
    if(!seq_char)
	    return printf("error malloc\n"), 1;
    seq_char_size = load_from_file("seq_test_intrinsic.txt", seq_char);

    seq_long = convert_to_binary(seq_char, seq_char_size);

    double cycles = 0.0;
    double boucles = 0.0;
    unsigned long long cnt = 0;

    // real detecting genes
    while(cycles <= 500)
    {
    	double before = (double)rdtsc();
    	cnt = test1(seq_char_size, seq_long);
    	double after = (double)rdtsc();

    	cycles += (after - before);
    	boucles++;
    }
    printf("real detecting genes : %lf cycles, %0.1lf boucles, genes cnt = %llu\n", cycles / boucles, boucles, cnt);

    // intrinsic detecting genes
    cycles = 0.0;
    boucles = 0.0;
    cnt = 0;
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        cnt = test2(seq_char_size, seq_long);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("intrinsic detecting genes : %lf cycles, %0.1lf boucles, genes cnt = %llu\n", cycles / boucles, boucles, cnt);

    free(seq_char);
    free(seq_long);

    return 0;
}
