#pragma once 

#define MAX_GENES  1024
// Number of bits in an integer
#define int_SIZE 63

typedef struct gene_map_s {

    //
    unsigned long long genes_counter;

    //Gene start position (AUG)
    unsigned long long* gene_start;

    //Gene stop position (UAA, UAG, UGA)
    unsigned long long* gene_end;

}gene_map_t;

typedef struct codon {
    char* codon;
    char* short_name;
    char* full_name;
    char symbol;
}codon;

typedef struct mutation_map {
    unsigned long* size;
    unsigned long *start_mut;
    unsigned long *end_mut;
}mutation_map;


/********** BINARIES FUNCTION **********/

int get_binary_value(const long int *seq_bin, const int pos);
long int* change_binary_value(long int *seq_bin, const int pos, const int value);
long int* set_binary_array(const char *array, const unsigned size);
long int* xor_binary_array(const long int *seq1, const unsigned array_size1,
                                const long int *seq2, const unsigned array_size2);
int popcount_binary_array(const long int *seq, const long int size);
long int mask_binary_array(const long int seq_bin, const long int pos_start, const long int size);
long int* get_piece_binary_array(const long int* seq_bin,  const long int pos_start, const long int size);


/******** DNA & GENES FUNCTION *********/

long int* convert_to_binary(const char* dna_seq, const unsigned size);
char* binary_to_dna(long int* bin_dna_seq, const unsigned size);
char* generating_mRNA(const long int* gene_seq, const long start_pos,const long int seq_size);
void detecting_genes(const long int *gene, const long int gene_size,
                     gene_map_t* gene_map);
char* generating_amino_acid_chain(const long int *gene_seq,const long int start_pos, const long int seq_size);
void detecting_mutations(const long int *gene_seq,const long int start_pos, const long int size_sequence,
                         mutation_map mut_m);
float calculating_matching_score(const long int *seq1, long int start_pos1,const int seq_size1,
                                 const long int *seq2, long int start_pos2,const int seq_size2);


