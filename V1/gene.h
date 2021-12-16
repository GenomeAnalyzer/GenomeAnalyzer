#pragma once 
#define MAX_GENES  1024

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

//extern struct gene_map_s gene_map;
short* convert_to_binary(char* dna_seq, unsigned size);
char* generating_mRNA(const unsigned int gene_seq [], const unsigned int seq_size);
void detecting_genes(const unsigned int gene [], const unsigned int gene_size, gene_map_t* gene_map);
char* generating_amino_acid_chain(char* seq, int seq_size, int codons_count, codon codons []);
void detecting_mutations(const unsigned short gene_seq [], const unsigned long size_sequence, mutation_map mut_m);
float calculating_matching_score(int sequence_size, int seq1 [], int seq2 []);
int binary_size_count(int b);
int hamming(int seq1, int seq2);