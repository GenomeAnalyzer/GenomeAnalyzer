#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.h"

//A   : 00
//T/U : 11
//C   : 10
//G   : 01

// UAA 110000
// UAG 110001
// UGA 110100

// ATC 001110

//////////////// Generating mRNA
/**
 * generating_mRNA - Convert a binary DNA sequence to a string mRNA sequence
 * dna_seq      array of binary
 * seq_len      size of dna_seq
 *
 * rna_seq      array of char
 */
char* generating_mRNA(int* dna_seq, int seq_len) {
    // Check the input argument
    if (!dna_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Create and check the output
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_len / 2));
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;
    // Parse the binary DNA sequence two by two
    for (int i = 0; i < seq_len; i += 2) {
        switch (dna_seq[i]) {
        case 0:
            if (dna_seq[i + 1] == 0)
                rna_seq[j] = 'A';
            else if (dna_seq[i + 1] == 1)
                rna_seq[j] = 'G';
            break;
        case 1:
            if (dna_seq[i + 1] == 0)
                rna_seq[j] = 'C';
            else if (dna_seq[i + 1] == 1)
                rna_seq[j] = 'U';
            break;
        default:
            return printf("ERROR: generating_mRNA: invalid value in DNA sequence\n"), NULL;
        }
        j++;
    }
    return rna_seq;
}

//////////////// Algo naif
void detecting_genes(unsigned int gene [], gene_map_t* gene_map) {
    //struct gene_map_s gene_map;
    gene_map->genes_counter = 0;

    int start_pos = -1;
    int stop_pos = -1;

    int i = 0;

    while ((i + 6) < 1000)     {


        if (start_pos == -1 && stop_pos == -1)     {

            // if (!(gene[i%32] & ( 1 << (i%32) ))
            if (gene[i] == 0 && gene[i + 1] == 0 && gene[i + 2] == 1 && gene[i + 3] == 1 && gene[i + 4] == 1 && gene[i + 5] == 0)         {
                start_pos = i;
            }
        }
        else
            if (start_pos != -1 && stop_pos == -1)     {
                if (((gene[i] == 1 && gene[i + 1] == 1 && gene[i + 2] == 0) && (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 0)) || (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 1) || (gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 0))          {
                    stop_pos = i;
                }
            }
            else {

                // if (start_pos != -1 && stop_pos != -1 )

                gene_map->gene_start[gene_map->genes_counter] = start_pos;
                gene_map->gene_end[gene_map->genes_counter] = stop_pos;

                gene_map->genes_counter++;

                start_pos = -1;
                stop_pos = -1;
            }

        i++;
    }

    // Algorithm:
        // ----------
        //     start_pos = lookup AUG in seq;
        // stop_pos = lookup(UAA or UAG or UGA) in seq + start_pos;
        // if (start_pos && stop_pos)   {
        //     gene_map.gene_start[gene_map.genes_counter) = start_pos;
        //     gene_map.gene_stop[gene_map.genes_counter) = stop_pos;
}

//////////////// Generating an amino acid chain (protein) 
void generating_amino_acid_chain() {

}


//////////////// Detecting probable mutation zones
/*
* in : gene_seq : sequence of the gene
* in : size_sequence : size of the sequence
* out : boolean
* The algorithm runs through a gene sequence and detects if there is a sequence with a high GC frequency (at least
* 1/5th of the gene sequence's size) it returns true, else it returns false.
* Precondition: gene_seq is of size size_sequence.
*/
bool detecting_mutations(const unsigned int gene_seq [], const unsigned long long size_sequence) {

    unsigned long long detect_mut = 0;  //Counting size of GC sequence

    //Read the sequence
    for (unsigned long long i = 0; i < size_sequence; i += 2) {

        //Increment detect_mut if find a C or G gene
        if (((__builtin_popcount(gene_seq[i]) == 0) && (__builtin_popcount(gene_seq[i + 1]) == 1)) ||
            ((__builtin_popcount(gene_seq[i]) == 1) && (__builtin_popcount(gene_seq[i + 1]) == 0))) {
            detect_mut++;
        }
        //Put detect_mut to 0 if find a A or T gene
        else {
            //Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= ((size_sequence / 2) / 5)) {
                return true;
            }
            detect_mut = 0;
        }
    }
    //Check if ending sequence is a probable mutation zone
    if (detect_mut >= ((size_sequence / 2) / 5)) {
        return true;
    }
    return false;
}


//////////////// Calculating the matching score of two sequences
void calculating_matching_score() {

}

//////////////// Hamming calculation
void hamming(int b1, int b2) {
    int x = b1 ^ b2;
    return __builtin_popcount(x);
}