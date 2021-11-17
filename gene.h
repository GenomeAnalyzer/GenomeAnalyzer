#pragma once 
#define MAX_GENES  1024

typedef struct gene_map_s{

   //
   unsigned long long genes_counter;

   //Gene start position (AUG)
   unsigned long long *gene_start;

   //Gene stop position (UAA, UAG, UGA)
   unsigned long long *gene_end;

}gene_map_t;

typedef struct codon{
    int number_of_codons;
    char *codons[3];
    char *short_name;
    char symbol;
    char *full_name;
}codon;

//extern struct gene_map_s gene_map;

void  detecting_genes(unsigned int gene[],gene_map_t *gene_map);
bool detecting_mutations(const unsigned int gene_seq[], const unsigned long long size_sequence);