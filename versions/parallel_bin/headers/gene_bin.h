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
#define REGSIZE 512
#define STEPCHAR 256
#else
#define REGSIZE 256
#define STEPCHAR 128
#endif

#include <mpi.h>
#include <sys/types.h>

typedef struct gene_map_s
{

    //
    unsigned long long genes_counter;

    // Gene start position (AUG)
    unsigned long long *gene_start;

    // Gene stop position (UAA, UAG, UGA)
    unsigned long long *gene_end;

} gene_map_t;

typedef struct mutation_map
{
    unsigned long *size;
    unsigned long *start_mut;
    unsigned long *end_mut;
} mutation_map;


typedef struct node
{
    long int *seq;
    int size;
    struct node *next;
    struct node *prev;

} node_t;

/***************************************/

// Initialisation of lookup table
typedef int lookuptable[2];

// Bit values according to ASCII code of nucleotides - 65
static lookuptable L[25] = {
    {0, 0},
    {1, 0},
    {1, 0},
    {0, 0},
    {-1, -1},
    {-1, -1},
    {0, 1},
    {0, 1},
    {-1, -1},
    {-1, -1},
    {0, 1},
    {-1, -1},
    {0, 0},
    {0, 0},
    {-1, -1},
    {-1, -1},
    {-1, -1},
    {0, 0},
    {1, 0},
    {1, 1},
    {-1, -1},
    {0, 0},
    {0, 0},
    {-1, -1},
    {1, 0}};

//******************************
// Code ASCII de A,G,C,T
static int bitstocharDNA[4] = {65, 71, 67, 84};
//******************************
// Code ASCII de A,G,C,U
static int bitstocharmRNA[4] = {65, 71, 67, 85};
//******************************
// Lookup Table Initialization
static char LUT[64] = {'K', 'K', 'N', 'N', 'R', 'R', 'S', 'S', 'T', 'T',
                       'T', 'T', 'I', 'M', 'I', 'I', 'E', 'E', 'D', 'D',
                       'G', 'G', 'G', 'G', 'A', 'A', 'A', 'A', 'V', 'V',
                       'V', 'V', 'Q', 'Q', 'H', 'H', 'R', 'R', 'R', 'R',
                       'P', 'P', 'P', 'P', 'L', 'L', 'L', 'L', 'O', 'O',
                       'Y', 'Y', 'O', 'W', 'C', 'C', 'S', 'S', 'S', 'S',
                       'L', 'L', 'F', 'F'};

/********** BINARIES FUNCTION **********/

int get_binary_value(const long int *seq_bin, const int pos);
long int *change_binary_value(long int *seq_bin, const int pos, const int value);
long int *set_binary_array(const char *array, const size_t size);
long int *xor_binary_array(const long int *seq1, const int array_size1,
                           const long int *seq2, const int array_size2);
int popcount_binary_array(const long int *seq, const long int size);
long int *get_piece_binary_array(const long int *seq_bin, const unsigned long long pos_start, const unsigned long long pos_stop);

/******** DNA & GENES FUNCTION *********/

long* convert_to_binary(const char *dna_seq, size_t size);
char *binary_to_dna(long int *bin_dna_seq, const unsigned size);
char *generating_mRNA(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos);
void detecting_genes(const long int *gene, const long int gene_size,
                     gene_map_t *gene_map);
char *generating_amino_acid_chain(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos);
void detecting_mutations(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos,
                         mutation_map mut_m);
float calculating_matching_score(long int *seq1, const int seq_size1,
                                 long int *seq2, const int seq_size2);

void launch();
