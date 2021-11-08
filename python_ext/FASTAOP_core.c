#include <stdio.h>
#include <stdlib.h>

int hamming(int b1, int b2){
    int x = b1 ^ b2;
    return __builtin_popcount(x);
}

char* RNA_generate(int seq_len, char* seq) {
    // for (unsigned long long i = 0; i < seq_len; i++)
    //     if (seq[i] == 'T')
    //         seq[i] = 'U';
    return seq;
}

//////////////// Detecting genes
int genes_detect(int var) {
    // struct gene_map_s {
    //     //
    //     unsigned long long genes_counter;
    //     //Gene start position (AUG)
    //     unsigned long long gene_start[MAX_GENES];
    //     //Gene stop position (UAA, UAG, UGA)
    //     unsigned long long gene_end[MAX_GENES];
    // };
    // gene_map_s gene_map;


// Algorithm:
    // ----------
    //     start_pos = lookup AUG in seq;
    // stop_pos = lookup(UAA or UAG or UGA) in seq + start_pos;
    // if (start_pos && stop_pos)   {
    //     gene_map.gene_start[gene_map.genes_counter) = start_pos;
    //     gene_map.gene_stop[gene_map.genes_counter) = stop_pos;
    return var;
}

//////////////// Generating an amino acid chain (protein) 
int amino_acid_chain_generate(int var) {
return var;
}

//////////////// Detecting probable mutation zones
int mutations_detect(int var) {
return var;
}

//////////////// Calculating the matching score of two sequences
int matching_score_calculate(int a, int b) {
    // return matching score
    return a+b;
}