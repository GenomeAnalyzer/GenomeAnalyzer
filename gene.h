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
    short codon[6];
    char* short_name;
    char* full_name;
    char symbol;
}codon;

//extern struct gene_map_s gene_map;

char* generating_mRNA(const unsigned int gene_seq [], const unsigned int seq_size);
void detecting_genes(const unsigned int gene [], const unsigned int gene_size, gene_map_t* gene_map);
// char* generating_amino_acid_chain(int seq_size, char* seq, const unsigned int gene_size, gene_map_t* gene_map, int codons_count, codon codons []);
// char* generating_amino_acid_chain(int seq_size, char* seq, int codons_count, codon codons []);
// short* compress(int size, char* dna_seq);
// char* generating_amino_acid_chain(int seq_size, char* seq);
void generating_amino_acid_chain(char** proteines, short seq_size, short* seq, gene_map_t gene_map, short codons_count, codon codons[]);
char codon_binary_to_aa(short* cdn, short codons_count, codon* codons);
bool detecting_mutations(const unsigned int gene_seq [], const unsigned long long size_sequence);
float calculating_matching_score(int sequence_size, int seq1 [], int seq2 []);
int binary_size_count(int b);
int hamming(int seq1, int seq2);