#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../headers/gene.h"

/**
 * in : dna_seq : array of char
 * in : size : size of dna_seq
 * out : seq : array of binary
 * Convert a char DNA sequence to its binary sequence
 */
unsigned short* convert_to_binary(char* dna_seq, unsigned size){
    unsigned temp = 0;

    //Create and check the output
    unsigned short *seq = calloc(sizeof(unsigned short),size); 
    if(!seq)
        return printf("ERROR: convert_to_binary: cannot allocate memory.\n"), NULL;

    //Parse the sequence, two nucleotides by two
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

            case 'R': // R = A or G, here R = A
                seq[temp] = 0;
                seq[temp+1] = 0; 
                temp += 2;
                break;

            case 'Y': // Y = C or T, here Y = C
                seq[temp] = 1;
                seq[temp+1] = 0; 
                temp += 2;
                break;

            case 'K': // K = G or T, here K = G
                seq[temp] = 0;
                seq[temp+1] = 1; 
                temp += 2;
                break;

            case 'M': // M = A or C, here M = A
            seq[temp] = 0;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            case 'S': // S = C or G, here S = C
            seq[temp] = 1;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            case 'W': // W = T or A, here W = T
            seq[temp] = 1;
            seq[temp+1] = 1; 
            temp += 2;
            break;

            case 'B': // B = not A, here B = C
            seq[temp] = 1;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            case 'D': // D = not C, here D = A
            seq[temp] = 0;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            case 'H': // H = not G, here H = A
            seq[temp] = 0;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            case 'V': // V = not T, here V = A
            seq[temp] = 0;
            seq[temp+1] = 0; 
            temp += 2;
            break;

            default:
                printf("Error: convert_to_binary: wrong letter in the sequence(%c).\nExit.\n",dna_seq[i]);
                seq = NULL;
                return seq;
        }
    }

    return seq;
}

/**
 * in : bin_dna_seq : unsigned short array - size must be 3*2 = 6.
 * out : aa : array of binary
 * Convert a binary aa sequence to its aa codon
 */
char* binary_to_aa(unsigned short* bin_dna_seq, unsigned size){
    //Check the input
    if (size % 2 != 0) {
        printf("Error: wrong binary size (%d). Must be odd.\nExit.\n",size);
        return NULL;
    }

    //Create and check the output
    char* codon = malloc(sizeof(char)*size/2);
    if (!codon)
        return printf("ERROR: binary_to_aa: cannot allocate memory.\n"), NULL;

    //Parse the array, two by two
    for (unsigned i = 0;i < size/2; i ++){
        if (bin_dna_seq[2*i] == 0 && bin_dna_seq[2*i+1] == 0)
            codon[i] = 'A';
        else if (bin_dna_seq[2*i] == 1 && bin_dna_seq[2*i+1] == 1)
            codon[i] = 'T';
        else if (bin_dna_seq[2*i] == 1 && bin_dna_seq[2*i+1] == 0)
            codon[i] = 'C';
        else if (bin_dna_seq[2*i] == 0 && bin_dna_seq[2*i+1] == 1)
            codon[i] = 'G';
    }

    return codon;
}

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
    rna_seq = malloc(sizeof(*rna_seq) * (seq_size / 2) + 2);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;
    // Parse the binary DNA sequence two by two
    for (unsigned i = 0; i < seq_size; i += 2) {
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
    rna_seq[j] = '\0';
    return rna_seq;
}

/*
* in : gene : sequence of genes
* in : gene_map : struct to map the genes
* out : void
* Detect if a gene exists in the sequence and insert it in the structure
*/
void detecting_genes(const unsigned long gene [], const unsigned int gene_size, gene_map_t* gene_map) {
    // Initialise the gene counter
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

    bool find = false;
    unsigned long long start_pos = -1;
    unsigned long long i = 0;

    //Parse the array
    while ((i + 6) <= gene_size) {   

            //If a 'start pos' and a 'stop pos' doesn't exist, search for AUG
            if (gene[i] == 0 && gene[i + 1] == 0 && gene[i + 2] == 1 
                && gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 1) {
            //if exist, it's the start of a gene
                start_pos = i;
                i += 6;
                find = true;
            }
            else{

            if (find ) {
                //if a 'start pos' exists , search for UAA / UAG / UGA
                if ((gene[i] == 1 && gene[i + 1] == 1 && gene[i + 2] == 0) 
                    && ((gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 0)
                        || (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 1)
                        || (gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 0))) {
                   //If exist, it's the end of a gene          
                   //If a 'start pos' and an 'stop pos' has been found,
                   //a gene exists so we save it in the struc
                    gene_map->gene_start[gene_map->genes_counter] = start_pos;
                    gene_map->gene_end[gene_map->genes_counter] = i+5;

                    gene_map->genes_counter++;

                    start_pos = -1;
                    i += 6;
                    find = false;
                }
                else
                    i += 2;
            }
            else
                i += 2;
            }
        
    }
}

/*
 * in : seq : original mRNA sequence.
 * out : char* : protein in symbols
 * The program parses the mRNA sequence, verify its length and if the first codon is a START codon.
*/
char* generating_amino_acid_chain(const unsigned short gene_seq [], const unsigned int seq_size) {
    // 1 codon = 3 nucleotides = 6 bits
    short codon_size = 6;

    // Check the input arguments
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;
    if(seq_size % 3 != 0)
        return NULL;

    // Create and check the output
    char* aa_seq = NULL;
    aa_seq = malloc(sizeof(*aa_seq) * (seq_size / codon_size) +1 );
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned temp = 0;

    for (unsigned int i = 0; i < seq_size; i += codon_size) {
        // The hash functions, takes the 6 bits, and transform the array into an integer.
        // The integer first char is a 2, for hash generation purposes.
        int hash = 2;
        for(unsigned k = i; k<i+codon_size; k++){
            hash = 10 * hash + gene_seq[k];
        }

        // Switch over the hash.
        switch(hash){
        case 2000000 :
            aa_seq[temp] = 'K';
            break;
        case 2000001 :
            aa_seq[temp] = 'K';
            break;
        case 2000010 :
            aa_seq[temp] = 'N';
            break;
        case 2000011 :
            aa_seq[temp] = 'N';
            break;
        case 2000100 :
            aa_seq[temp] = 'R';
            break;
        case 2000101 :
            aa_seq[temp] = 'R';
            break;
        case 2000110 :
            aa_seq[temp] = 'S';
            break;
        case 2000111 :
            aa_seq[temp] = 'S';
            break;
        case 2001000 :
            aa_seq[temp] = 'T';
            break;
        case 2001001 :
            aa_seq[temp] = 'T';
            break;
        case 2001010 :
            aa_seq[temp] = 'T';
            break;
        case 2001011 :
            aa_seq[temp] = 'T';
            break;
        case 2001100 :
            aa_seq[temp] = 'I';
            break;
        case 2001101 :
            aa_seq[temp] = 'M';
            break;
        case 2001110 :
            aa_seq[temp] = 'I';
            break;
        case 2001111 :
            aa_seq[temp] = 'I';
            break;
        case 2010000 :
            aa_seq[temp] = 'E';
            break;
        case 2010001 :
            aa_seq[temp] = 'E';
            break;
        case 2010010 :
            aa_seq[temp] = 'D';
            break;
        case 2010011 :
            aa_seq[temp] = 'D';
            break;
        case 2010100 :
            aa_seq[temp] = 'G';
            break;
        case 2010101 :
            aa_seq[temp] = 'G';
            break;
        case 2010110 :
            aa_seq[temp] = 'G';
            break;
        case 2010111 :
            aa_seq[temp] = 'G';
            break;
        case 2011000 :
            aa_seq[temp] = 'A';
            break;
        case 2011001 :
            aa_seq[temp] = 'A';
            break;
        case 2011010 :
            aa_seq[temp] = 'A';
            break;
        case 2011011 :
            aa_seq[temp] = 'A';
            break;
        case 2011100 :
            aa_seq[temp] = 'V';
            break;
        case 2011101 :
            aa_seq[temp] = 'V';
            break;
        case 2011110 :
            aa_seq[temp] = 'V';
            break;
        case 2011111 :
            aa_seq[temp] = 'V';
            break;
        case 2100000 :
            aa_seq[temp] = 'Q';
            break;
        case 2100001 :
            aa_seq[temp] = 'Q';
            break;
        case 2100010 :
            aa_seq[temp] = 'H';
            break;
        case 2100011 :
            aa_seq[temp] = 'H';
            break;
        case 2100100 :
            aa_seq[temp] = 'R';
            break;
        case 2100101 :
            aa_seq[temp] = 'R';
            break;
        case 2100110 :
            aa_seq[temp] = 'R';
            break;
        case 2100111 :
            aa_seq[temp] = 'R';
            break;
        case 2101000 :
            aa_seq[temp] = 'P';
            break;
        case 2101001 :
            aa_seq[temp] = 'P';
            break;
        case 2101010 :
            aa_seq[temp] = 'P';
            break;
        case 2101011 :
            aa_seq[temp] = 'P';
            break;
        case 2101100 :
            aa_seq[temp] = 'L';
            break;
        case 2101101 :
            aa_seq[temp] = 'L';
            break;
        case 2101110 :
            aa_seq[temp] = 'L';
            break;
        case 2101111 :
            aa_seq[temp] = 'L';
            break;
        case 2110000 :
            aa_seq[temp] = 'O';
            break;
        case 2110001 :
            aa_seq[temp] = 'O';
            break;
        case 2110010 :
            aa_seq[temp] = 'Y';
            break;
        case 2110011 :
            aa_seq[temp] = 'Y';
            break;
        case 2110100 :
            aa_seq[temp] = 'O';
            break;
        case 2110101 :
            aa_seq[temp] = 'W';
            break;
        case 2110110 :
            aa_seq[temp] = 'C';
            break;
        case 2110111 :
            aa_seq[temp] = 'C';
            break;
        case 2111000 :
            aa_seq[temp] = 'S';
            break;
        case 2111001 :
            aa_seq[temp] = 'S';
            break;
        case 2111010 :
            aa_seq[temp] = 'S';
            break;
        case 2111011 :
            aa_seq[temp] = 'S';
            break;
        case 2111100 :
            aa_seq[temp] = 'L';
            break;
        case 2111101 :
            aa_seq[temp] = 'L';
            break;
        case 2111110 :
            aa_seq[temp] = 'F';
            break;
        case 2111111 :
            aa_seq[temp] = 'F';
            break;

        default:
            return printf("ERROR: generating_amino_acid_chain: invalid value (%d) in RNA sequence\n", hash), NULL;
        }
    
        temp++;
    }
    aa_seq[temp] = '\0';
    return aa_seq;
}

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
    unsigned long long tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones


    //If problem in allocation
    if(!mut_m.start_mut || !mut_m.end_mut || !mut_m.size){
        mut_m.start_mut = malloc(5 * sizeof(mut_m.start_mut));
        mut_m.start_mut = malloc(5 * sizeof(mut_m.end_mut));
        mut_m.start_mut = malloc(5 * sizeof(mut_m.size));
    }

    //Parse the sequence
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

    //Find the greater sequence and rename both
    const unsigned short* s1 = sequence_size1 >= sequence_size2 ? seq1 : seq2;
    const unsigned short* s2 = sequence_size1 >= sequence_size2 ? seq2 : seq1;
    int diff_size = sequence_size1 >= sequence_size2 ? sequence_size1 - sequence_size2 : sequence_size2 - sequence_size1;
    int max_size = sequence_size1 >= sequence_size2 ? sequence_size1 : sequence_size2;
    int count = 0;

    // For the greater part of the greater sequence, count the number of '1'
    for (int i = 0; i < diff_size; i++)
        if (s1[i]) ++count;

    // Count the number of '1', after aplying the xor operation
    for (int i = diff_size; i < max_size; ++i) {
        if (s1[i] ^ s2[i - diff_size]) ++count;
    }

    // Compute and return the percentage
    float y = ((float)count * 100.0) / (float)max_size;
    return 100.0 - y;
}
