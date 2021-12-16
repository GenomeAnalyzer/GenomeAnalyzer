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

//////////////// Convert to binary
/**
 * in : dna_seq : array of char
 * in : size : size of dna_seq
 * out : seq : array of binary
 * Convert a char DNA sequence to its binary sequence
 */
unsigned short* convert_to_binary(char* dna_seq, unsigned size){
    unsigned i = 0;
    unsigned temp = 0;

    short *seq = malloc(sizeof(short)*size); 

    for (unsigned i = 0;i < size/2; i ++){
        switch(dna_seq[i]){
            case 'A':
                seq[temp] = 0;
                seq[temp+1] = 0; 
                temp += 2;

                break;
            
            case 'T':
                seq[temp] = 1;
                seq[temp+1] = 1;
                temp += 2;

                break;
            
            case 'C':
                seq[temp] = 1;
                seq[temp+1] = 0;
                temp += 2;

                break;
            case 'G':
                seq[temp] = 0;
                seq[temp+1] = 1;
                temp += 2;

                break;
            case 'N':
                seq[temp] = 0;
                seq[temp+1] = 0; 
                temp += 2;

                break;
            default:
                printf("Error: wrong letter in the sequence(%c).\nExit.\n",dna_seq[i]);

                return seq;
        }
    }

    return seq;
}

//////////////// Convert binary aa to codon
/**
 * in : bin_dna_seq : unsigned short array - size must be 3*2 = 6.
 * out : aa : array of binary
 * Convert a binary aa sequence to its aa codon
 */
char* binary_to_aa(unsigned short* bin_dna_seq, unsigned size){
    if (size%2 != 0) {
        printf("Error: wrong binary size (%d). Must be odd.\nExit.\n",size);
        return NULL;
    }

    // char *codon[size/2];
    char* codon = malloc(sizeof(char)*size/2);
    for (unsigned i = 0;i < size/2; i ++){
        if (bin_dna_seq[2*i] == 0 && bin_dna_seq[2*i+1] == 0)
            codon[i] = 'A'; // Assuming "N" is "A", as it's done in convert_to_binary
        else if (bin_dna_seq[2*i] == 1 && bin_dna_seq[2*i+1] == 1)
            codon[i] = 'T';
        else if (bin_dna_seq[2*i] == 1 && bin_dna_seq[2*i+1] == 0)
            codon[i] = 'C';
        else if (bin_dna_seq[2*i] == 0 && bin_dna_seq[2*i+1] == 1)
            codon[i] = 'G';
    }

    return codon;
}

//////////////// Generating mRNA
/**
 * in : gene_seq : array of binary
 * in : seq_len : size of dna_seq
 * out : rna_seq : array of char
 * Convert a binary DNA sequence to a string mRNA sequence
 */
char* generating_mRNA(const unsigned short gene_seq [], const unsigned int seq_size) {
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Create and check the output
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_size / 2) + 1);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;
    // Parse the binary DNA sequence two by two
    for (int i = 0; i < seq_size; i += 2) {
        switch (gene_seq[i]) {
        case 0:
            if (gene_seq[i + 1] == 0)
                // A = 00
                rna_seq[j] = 'A';
            else if (gene_seq[i + 1] == 1)
                // G = 01
                rna_seq[j] = 'G';
            break;
        case 1:
            if (gene_seq[i + 1] == 0)
                // C = 10
                rna_seq[j] = 'C';
            else if (gene_seq[i + 1] == 1)
                // U = 11
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
/*
* in : gene : sequence of genes
* in : gene_map : struct to map the genes
* out : void
* Detect if a gene exists in the sequence and insert it in the structure
*/
void detecting_genes(const unsigned int gene [], const unsigned int gene_size, gene_map_t* gene_map) {
    //struct gene_map_s gene_map;
    gene_map->genes_counter = 0;

    // Check if memory ever have been allocated, do it if not
    if(!gene_map->gene_start || !gene_map->gene_end){
        gene_map->gene_start = malloc(sizeof(*gene_map->gene_start) * MAX_GENES);
        gene_map->gene_end = malloc(sizeof(*gene_map->gene_end) * MAX_GENES);

        if (!gene_map->gene_start || !gene_map->gene_end){
            printf("ERROR: detecting_genes: cannot allocate memory\n");
            return;
        }
    }


    int start_pos = -1;

    int i = 0;

    while ((i + 6) <= gene_size) {   

            //If a start pos and a stop pos doesn't exist, search for AUG
            // if (!(gene[i%32] & ( 1 << (i%32) ))
            if (gene[i] == 0 && gene[i + 1] == 0 && gene[i + 2] == 1 
                && gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 1) {
            //if atc, it's the start of a gene
                start_pos = i;
                i += 6;
            }
            else{

            if (start_pos != -1 ) {
                //if a start pos exists , search for UAA / UAG / UGA
                if ((gene[i] == 1 && gene[i + 1] == 1 && gene[i + 2] == 0) 
                    && ((gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 0)
                        || (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 1)
                        || (gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 0))) {
                   //It's the end of a gene          
                   //If a start pos and an stop pos has been found, a gene exists so we save it in the struc
                    gene_map->gene_start[gene_map->genes_counter] = start_pos;
                    gene_map->gene_end[gene_map->genes_counter] = i;

                    gene_map->genes_counter++;

                    start_pos = -1;
                    i += 6;
                }
                else
                    i += 2;
            }
            else
                i += 2;
            }
        
    }
}
//////////////// Generating an amino acid chain (protein) 
/*
 * in : seq : original mRNA sequence.
 * out : char* : protein in symbols
 * The program parses the mRNA sequence, verify its length and if the first codon is a START codon.
*/
char* generating_amino_acid_chain(const unsigned short gene_seq [], const unsigned int seq_size) {
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;

    // Create and check the output
    char* aa_seq = NULL;
    aa_seq = malloc(sizeof(*aa_seq) * (seq_size / 6) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned int j = 0;

    // Parse the binary DNA sequence six by six
   for (unsigned int i = 0; i < seq_size; i += 6) {

        // If Axx
        if (gene_seq[i] == 0 &&  gene_seq[i + 1] == 0){

            // If AAx
            if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 0) {

                if (gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // AAA = K
                    aa_seq[j] = 'K';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // AAC = N
                    aa_seq[j] = 'N';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // AAG = K
                    aa_seq[j] = 'K';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // AAU = N
                    aa_seq[j] = 'N';
            }

            // If ACx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 0) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // ACA = T
                    aa_seq[j] = 'T';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // ACC = T
                    aa_seq[j] = 'T';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // ACG = T
                    aa_seq[j] = 'T';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // ACU = T
                    aa_seq[j] = 'T';
            }

            // If AGx
            else if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // AGA = R
                    aa_seq[j] = 'R';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // AGC = S
                    aa_seq[j] = 'S';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // AGG = R
                    aa_seq[j] = 'R';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // AGU = S
                    aa_seq[j] = 'S';
            }

            // If AUx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // AUA = I
                    aa_seq[j] = 'I';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // AUC = I
                    aa_seq[j] = 'I';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // AUG = Start = M
                    aa_seq[j] = 'M';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // AUU = I
                    aa_seq[j] = 'I';
            }
        }

        // If Cxx
        else if (gene_seq[i] == 1 &&  gene_seq[i + 1] == 0){

            // If CAx
            if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 0) {

                if (gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // CAA = Q
                    aa_seq[j] = 'Q';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // CAC = H
                    aa_seq[j] = 'H';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // CAG = Q
                    aa_seq[j] = 'Q';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // CAU = N
                    aa_seq[j] = 'H';
            }

            // If CCx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 0) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // CCA = P
                    aa_seq[j] = 'P';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // CCC = T
                    aa_seq[j] = 'P';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // CCG = T
                    aa_seq[j] = 'P';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // CCU = T
                    aa_seq[j] = 'P';
            }

            // If CGx
            else if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // CGA = R
                    aa_seq[j] = 'R';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // CGC = R
                    aa_seq[j] = 'R';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // CGG = R
                    aa_seq[j] = 'R';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // CGU = R
                    aa_seq[j] = 'R';
            }

            // If CUx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // CUA = L
                    aa_seq[j] = 'L';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // CUC = L
                    aa_seq[j] = 'L';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // CUG = I
                    aa_seq[j] = 'L';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // CUU = L
                    aa_seq[j] = 'L';
            }
        }

        // If Gxx
        else if (gene_seq[i] == 0 &&  gene_seq[i + 1] == 1){

            // If GAx
            if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 0) {

                if (gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // GAA = E
                    aa_seq[j] = 'E';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // GAC = D
                    aa_seq[j] = 'D';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // GAG = E
                    aa_seq[j] = 'E';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // GAU = D
                    aa_seq[j] = 'D';
            }

            // If GCx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 0) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // GCA = A
                    aa_seq[j] = 'A';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // GCC = T
                    aa_seq[j] = 'A';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // GCG = T
                    aa_seq[j] = 'A';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // GCU = T
                    aa_seq[j] = 'A';
            }

            // If GGx
            else if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // CGA = G
                    aa_seq[j] = 'G';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // GGC = G
                    aa_seq[j] = 'G';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // GGG = G
                    aa_seq[j] = 'G';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // GGU = R
                    aa_seq[j] = 'G';
            }

            // If GUx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // GUA = V
                    aa_seq[j] = 'V';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // GUC = V
                    aa_seq[j] = 'V';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // GUG = V
                    aa_seq[j] = 'V';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // GUU = V
                    aa_seq[j] = 'V';
            }
        }

        // If Uxx
        else if (gene_seq[i] == 1 &&  gene_seq[i + 1] == 1){

            // If UAx
            if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 0) {

                if (gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // UAA = Stop = O
                    aa_seq[j] = 'O';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // UAC = Y
                    aa_seq[j] = 'Y';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // UAG = Stop = O
                    aa_seq[j] = 'O';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // UAU = Y
                    aa_seq[j] = 'Y';
            }

            // If UCx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 0) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // UCA = S
                    aa_seq[j] = 'S';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // UCC = S
                    aa_seq[j] = 'S';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // UCG = S
                    aa_seq[j] = 'S';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // UCU = S
                    aa_seq[j] = 'S';
            }

            // If UGx
            else if (gene_seq[i + 2] == 0 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // UGA = Stop = O
                    aa_seq[j] = 'O';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // UGC = C
                    aa_seq[j] = 'C';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // UGG = W
                    aa_seq[j] = 'W';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // UGU = C
                    aa_seq[j] = 'C';
            }

            // If UUx
            else if (gene_seq[i + 2] == 1 &&  gene_seq[i + 3] == 1) {

                if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 0)
                    // UUA = L
                    aa_seq[j] = 'L';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 0)
                    // UUC = F
                    aa_seq[j] = 'F';

                else if(gene_seq[i + 4] == 0 &&  gene_seq[i + 5] == 1)
                    // UUG = L
                    aa_seq[j] = 'L';

                else if(gene_seq[i + 4] == 1 &&  gene_seq[i + 5] == 1)
                    // UUU = F
                    aa_seq[j] = 'F';
            }
        }
        else{
            return printf("ERROR: generating_amino_acid_chain: invalid value in RNA sequence\n"), NULL;
        }

        j++;
    }

    return aa_seq;
}

//////////////// Detecting probable mutation zones
/*
* in : gene_seq : sequence of the gene
* in : size_sequence : size of the sequence
* out : boolean
* The algorithm runs through a gene sequence and detects if there is a sequence with a high GC frequency
* (at least 1/5th of the gene sequence's size) it returns true, else it returns false.
* Precondition: gene_seq is of size size_sequence.
*/
void detecting_mutations(const unsigned short gene_seq [], const unsigned long size_sequence,
    mutation_map mut_m) {
    unsigned long detect_mut = 0;  //Counting size of GC sequence
    unsigned short tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones

    //Read the sequence
    for (unsigned long i = 0; i < size_sequence; i += 2) {
        //Increment detect_mut if find a C or G nucl
        if (((gene_seq[i] == 0) && (gene_seq[i + 1] == 1)) ||
            ((gene_seq[i] == 1) && (gene_seq[i + 1]) == 0)) {
            if(detect_mut == 0){tmp_start_mut = i;}
            detect_mut+=2;
        }
        //Put detect_mut to 0 if find a A or T nucl
        else {
            //Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= (size_sequence / 5)) {
                mut_m.start_mut[cmp] = tmp_start_mut;
                mut_m.end_mut[cmp] = i-1;
                mut_m.size[cmp] = detect_mut-1;
                cmp++;
            }
            detect_mut = 0;
        }
    }
    //Check if ending sequence is a probable mutation zone
    if (detect_mut >= (size_sequence / 5)) {
        mut_m.start_mut[cmp] = tmp_start_mut;
        mut_m.end_mut[cmp] = size_sequence-1;
        mut_m.size[cmp] = detect_mut-1;
    }
}


//////////////// Calculating the matching score of two sequences
/*
 * in : seq1 : first sequence in binary
 * in : seq2 : second sequence in binary
 * out : float : matching score in %
 * The algorithms runs the hamming distance between two binary sequences, and return their matching score in %
*/
float calculating_matching_score(const unsigned short seq1 [], const int sequence_size1,
                                 const unsigned short seq2 [], const int sequence_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    // If the sequences don't have the same size, do (with x = 1 or 0):

    //  xxxxxxxxx
    // ^
    //  000xxxxxx
    // -----------
    //  xxxxxxxxx

    // And, 0 ^ x = x

    const unsigned short* s1 = sequence_size1 >= sequence_size2 ? seq1 : seq2;
    const unsigned short* s2 = sequence_size1 >= sequence_size2 ? seq2 : seq1;
    int diff_size = sequence_size1 >= sequence_size2 ? sequence_size1 - sequence_size2 : sequence_size2 - sequence_size1;
    int max_size = sequence_size1 >= sequence_size2 ? sequence_size1 : sequence_size2;
    int count = 0;

    for (int i = 0; i < diff_size; i++)
        if (s1[i]) ++count;


    for (int i = diff_size; i < max_size; ++i) {
        if (s1[i] ^ s2[i - diff_size]) ++count;
    }

    float y = ((float)count * 100.0) / (float)max_size;
    return 100.0 - y;
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
 * in : seq1 : first sequence in binary
 * in : seq2 : second sequence in binary
 * out : int : count of bits
 * The algorithms calculate the hamming distance between two binary sequences
*/
int hamming(int seq1, int seq2) {
    int x = seq1 ^ seq2;
    return __builtin_popcount(x);
}