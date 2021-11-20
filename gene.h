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

//extern struct gene_map_s gene_map;

int* convert_to_binary(char* dna_seq);
char* generating_mRNA(int* dna_seq, int seq_len);
void  detecting_genes(unsigned int gene [], gene_map_t* gene_map);
char* generating_amino_acid_chain(char* seq, int seq_size, int codons_count, codon codons []);
bool detecting_mutations(const unsigned int gene_seq [], const unsigned long long size_sequence);
float calculating_matching_score(int sequence_size, int seq1 [], int seq2 []);
int binary_size_count(int b);
int hamming(int seq1, int seq2);