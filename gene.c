#include <stdio.h>
#include <stdlib.h>
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

    while ((i + 6) < 1000) {

        if (start_pos == -1 && stop_pos == -1) {

            // if (!(gene[i%32] & ( 1 << (i%32) ))
            if (gene[i] == 0 && gene[i + 1] == 0 && gene[i + 2] == 1 && gene[i + 3] == 1 && gene[i + 4] == 1 && gene[i + 5] == 0) {
                start_pos = i;
            }
        }
        else
            if (start_pos != -1 && stop_pos == -1) {
                if (((gene[i] == 1 && gene[i + 1] == 1 && gene[i + 2] == 0) && (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 0)) || (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 1) || (gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 0)) {
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
/*
 * in : seq : original mRNA sequence.
 * out : char* : protein in symbols
*/
// not working
char* generating_amino_acid_chain(char* seq) {
    // codon* codons = NULL;
    // codons = malloc(sizeof(codon) * 21);
    // if (!codons)
    //     return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;
    codon codons[21];

    char* filename = "codons.txt";

    // would have been much easier in python with the .json file ._.
    FILE* fp = NULL;
    fp = fopen(filename, "r");
    if (!fp)
        return printf("ERROR: generating_amino_acid_chain: cannot open %s\n", filename), NULL;

    char buffer[30];
    char buffer2[30];
    char buffer3[1];
    char buffer4[30];
    char* current_codon_name = NULL;
    int codon_number = 0;

    for (int i = 0; i < 21; i++) {
        codons[i].number_of_codons = 0;
        codons[i].codons[0] = NULL;
        codons[i].short_name = NULL;
        codons[i].symbol = NULL;
        codons[i].full_name = NULL;
    }


    char file_contents[50];
    int i = -1;

    // First line is header
    // for (int i = 0; i < 50; i++) {
        // fscanf(fp, "%s%s%s%s", &buffer, &buffer2, &buffer3, &buffer4);
    current_codon_name = "K";
    while (fscanf(fp, "%[^\n] ", file_contents) != EOF) {
        i++;
        if (i == 0) continue;

        printf("> %s\n", file_contents);
        // sleep(1);
        // if (i==0) continue;
        sscanf(file_contents, "%s%s%s%s", &buffer, &buffer2, &buffer3, &buffer4);
        printf("%s_%s_%s_%s\n", buffer, buffer2, buffer3, buffer4);

        // if (i == 1) current_codon_name = buffer3;

        printf("ahhh\n");
        codons[codon_number].codons[codons[codon_number].number_of_codons] = buffer;
        printf("ahhh\n");
        codons[codon_number].number_of_codons++;
        printf("ahhh\n");
        // printf("strcmp : %d _ %s _ %s\n\n", (buffer3 == current_codon_name), buffer3, current_codon_name);
        printf("ahhh\n");
        if (buffer3 != current_codon_name) {
            codons[codon_number].short_name = buffer2;
            codons[codon_number].symbol = buffer3;
            codons[codon_number].full_name = buffer4;

            // printf("Full_name : %s",codons[codon_number].full_name);
            printf("Full_name : %s\n\tNumber of codons : %d\n\tCodon : %s\n\tshort_name : %s\n\tSymbol : %c\t\n",
                codons[codon_number].full_name, codons[codon_number].number_of_codons, codons[codon_number].codons[0], codons[codon_number].short_name, codons[codon_number].symbol);
            // printf("Full_name : %s\n\tNumber of codons : %d\n\tCodon : %s\n\tshort_name : %s\n\tSymbol : %s\t\n", 
            // codons[codon_number].full_name, codons[codon_number].number_of_codons, codons[codon_number].codons[0], codons[codon_number].short_name, codons[codon_number].symbol);
            codon_number++;
        }
        // strcpy(current_codon_name, buffer3);
        current_codon_name = buffer3[0];
    }
    printf("at least en of\n");
    codons[codon_number].short_name = buffer2;
    codons[codon_number].symbol = buffer3;
    codons[codon_number].full_name = buffer4;

    printf("at least en of for\n");

    for (int i = 0; i < 21; i++) {
        printf("Full_name : %s\n\tNumber of codons : %d\n\tCodon : %s\n\tshort_name : %s\n\tSymbol : %s\t\n",
            codons[i].full_name, codons[i].number_of_codons, codons[i].codons[0], codons[i].short_name, codons[i].symbol);
    }

    fclose(fp);
    free(codons);
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
/*
 * in : seq1 : First sequence in binary
 * in : seq2 : Second sequence in binary
 * out : float : matching score in %
 * The algorithms runs the hamming distance between two binary sequences, and return their matching score in %
*/

float calculating_matching_score(int sequence_size, int seq1 [], int seq2 []) {
    // Check outside of function that seq1 and seq2 are the same size

    int total_size_sequence = 0;
    int total_hamming_distance = 0;

    for (int i = 0; i < sequence_size, i++){
        int size = binary_size_count(seq1[i]);
        // Make check outside of function ?
        if (size != binary_size_count(seq2[i])) {
            printf("ERROR: generating_mRNA: wrong size sequence\n");
            return -1.0;
        }
        total_size_sequence += size;
        total_hamming_distance += hamming(seq1, seq2);
    }

    return 100 * total_hamming_distance / total_size_sequence;
}

//////////////// Counting binary size
/*
 * in : b : binary number
 * out : int : count of bits
 * The algorithm return the count of bits the binary integer is.
*/
int binary_size_count(int b) {
    int size = 0;
    while (b) {
        size++;
        b >>= 1;
    }
    return size;
}

//////////////// Hamming calculation
/*
 * in : seq1 : First sequence in binary
 * in : seq2 : Second sequence in binary
 * out : int : count of bits
 * The algorithms calculate the hamming distance between two binary sequences
*/
int hamming(int seq1, int seq2) {
    int x = seq1 ^ seq2;
    return __builtin_popcount(x);
}