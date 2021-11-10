#include <stdio.h>
#include <stdlib.h>


//////////////// Generating mRNA
/**
 * generating_mRNA - Convert a binary DNA sequence to a string mRNA sequence
 * dna_seq      array of binary
 * seq_len      size of dna_seq
 * 
 * rna_seq      array of char
 */
char *generating_mRNA(int *dna_seq, int seq_len){
    // Check the input argument
    if(!dna_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Create and check the output
    char *rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_len / 2));
    if(!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;
    // Parse the binary DNA sequence two by two
    for (unsigned long long i = 0; i < seq_len; i += 2){
        switch(dna_seq[i]){
            case 0:
                if(dna_seq[i + 1] == 0)
                    rna_seq[j] = 'A';
                else if(dna_seq[i + 1] == 1)
                    rna_seq[j] = 'G';
                break;
            case 1:
                if(dna_seq[i + 1] == 0)
                    rna_seq[j] = 'C';
                else if(dna_seq[i + 1] == 1)
                    rna_seq[j] = 'U';
                break;
            default:
                return printf("ERROR: generating_mRNA: invalid value in DNA sequence\n"), NULL;
        }
        j++;
    }
    return rna_seq;
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
void detecting_mutations(){

}

//////////////// Calculating the matching score of two sequences
void calculating_matching_score(){

}

//////////////// Hamming calculation
void hamming(int b1, int b2){
    int x = b1 ^ b2;
    return __builtin_popcount(x);
}
