#include "gene.h"
#include <stdio.h>
#include <stdlib.h>

//////////////// Generating mRNA
void generating_mRNA(){
    // for (unsigned long long i = 0; i < seq_len; i++)
    //     if (seq[i] == 'T')
    //         seq[i] = 'U';
}


//////////////// Detecting genes
void detecting_genes(){
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
}

//////////////// Generating an amino acid chain (protein) 
void generating_amino_acid_chain(){

}

//////////////// Detecting probable mutation zones
void detecting_mutations(unsigned long long size_sequence, unsigned int gene_seq[], unsigned int* mut_zone){

    unsigned long long detect_mut = 0;  //Counting size of GC sequence
    unsigned long long start_mut = -1;  //First gene of GC sequence
    unsigned long long end_mut = -1;    //Last gene of GC sequence

    //
    for(unsigned long long i = 0; i < size_sequence ; i++){

        //Increment detect_mut if find a C or G gene
        if(gene_seq[i] == G || gene_seq[i] == C){
            //
            if(detect_mut == 0){
                start_mut = i;
            }
            end_mut = i;
            detect_mut++;
        }
        //Put detect_mut to 0 if find a A or T gene
        else{
            //GC sequence is a probable mutation zone
            if(detect_mut == (size_sequence / 5)){
                mut_zone[0] = start_mut;
                mut_zone[1] = end_mut;
                return 0;
            }
            detect_mut = 0;
        }
    }
    //Ending sequence is a probable mutation zone
    if(detect_mut == (size_sequence / 5)){
        mut_zone[0] = start_mut;
        mut_zone[1] = end_mut;
    }
    return 0;
}

//////////////// Calculating the matching score of two sequences
void calculating_matching_score(){

}

//////////////// Hamming calculation
void hamming(int b1, int b2){
    int x = b1 ^ b2;
    return __builtin_popcount(x);
}