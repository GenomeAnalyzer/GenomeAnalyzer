#pragma once 

#define MAX_GENES 1024
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

typedef struct mutation_map {
    unsigned long* size;
    unsigned long *start_mut;
    unsigned long *end_mut;
}mutation_map;





/******** DNA & GENES FUNCTION *********/

long int* convert_to_binary(const char* dna_seq, const unsigned size);
char* binary_to_dna(_Bool * bin_dna_seq, const unsigned size);
char* generating_mRNA(const _Bool* gene_seq,const long int seq_size);
void detecting_genes(const _Bool *gene, const long int gene_size,
                     gene_map_t* gene_map);
char* generating_amino_acid_chain(const _Bool *gene_seq, const long int seq_size);
void detecting_mutations(const _Bool *gene_seq, const long int size_sequence,
                         mutation_map mut_m);
float calculating_matching_score(const _Bool *seq1, const int seq_size1,
                                 const _Bool *seq2, const int seq_size2);


