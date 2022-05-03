#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <stdint.h>
#include <x86intrin.h>
#include <immintrin.h>

#define MAX_GENES 1024
// Number of bits in an integer
#define int_SIZE 63
#define OFFSET_TABLE 65

#ifdef __AVX512__
#define STEPCHAR 256
#else
#define STEPCHAR 128
#endif

typedef struct gene_map_s {

    //
    unsigned long long genes_counter;

    //Gene start position (AUG)
    unsigned long long* gene_start;

    //Gene stop position (UAA, UAG, UGA)
    unsigned long long* gene_end;

}gene_map_t;

typedef struct mutation_map {
    unsigned long* size;
    unsigned long *start_mut;
    unsigned long *end_mut;
}mutation_map;

typedef struct mm_array_s
{
#ifdef __AVX512__
    __m512i reg;
#else
    __m256i reg;
#endif
} mm_array_t;

//Initialisation of lookup table
    typedef int lookuptable[2];

    //Bit values according to ASCII code of nucleotides - 65
    static lookuptable L[25] = {
        {0,0},
        {1,0},
        {1,0},
        {0,0},
        {-1,-1},
        {-1,-1},
        {0,1},
        {0,1},
        {-1,-1},
        {-1,-1},
        {0,1},
        {-1,-1},
        {0,0},
        {0,0},
        {-1,-1},
        {-1,-1},
        {-1,-1},
        {0,0},
        {1,0},
        {1,1},
        {-1,-1},
        {0,0},
        {0,0},
        {-1,-1},
        {1,0}
    };


/********** BINARIES FUNCTION **********/

int get_binary_value(const long int *seq_bin, const int pos);


/******** DNA & GENES FUNCTION *********/

void my_convert_to_binary(mm_array_t *seq_bin, const uint64_t seq_bin_size, const char* seq_char, const uint64_t seq_char_size);
char* binary_to_dna(long int* bin_dna_seq, const unsigned size);
char* generating_mRNA(const long int* gene_seq, const long start_pos,const long int seq_size);
void detecting_genes(const long int *gene, const long int gene_size, gene_map_t* gene_map);
char* generating_amino_acid_chain(const long int *gene_seq,const long int start_pos, const long int seq_size);
void detecting_mutations(const long int *gene_seq,const long int start_pos, const long int size_sequence,
                         mutation_map mut_m);
float calculating_matching_score(const long int *seq1, long int start_pos1,const int seq_size1,
                                 const long int *seq2, long int start_pos2,const int seq_size2);


