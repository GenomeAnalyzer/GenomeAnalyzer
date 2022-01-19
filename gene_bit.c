#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "gene_bit.h"


/***************************************/
/********** BINARIES FUNCTION **********/
/***************************************/





/**
 * Xor two binary array sequences.
 * 
 * in : seq_bin1 : first sequence in binary array format to xor
 * in : seq_size1 : seq_bin1 length (number total of used bits)
 * in : seq_bin2 : first sequence in binary array format to xor
 * in : seq_size2 : seq_bin2 length (number total of used bits)
 * out : xor : binary array sequence resulting from the xor operation between seq1 and seq2
 * 
 * Iterates over sequences value per value. 
 * Xor the two sequences values since its value is the same length.
 * If one sequence is larger than the other, shift the last value of the smaller sequence to xor it with the other value.
 * Values from the largest binary array are assigned to the xor result. (x^0 = x)
 */
long int* xor_binary_array(_Bool* const seq_bin1, const unsigned seq_size1,
                           _Bool * const seq_bin2, const unsigned seq_size2){

    // size of the binary array type used
    long int intsize = int_SIZE + 1;

    long int* s1, * s2;
    long int ss1, ss2;

    // Find the greater binary array, and rename them
    if (seq_size1 >= seq_size2) { // if s1 is greater or equal than s2
        s1 = seq_bin1;
        s2 = seq_bin2;
        ss1 = seq_size1 ;
        ss2 = seq_size2 ;
    }
    else { // else if s2 is greater or equal than s1
        s1 = seq_bin2;
        s2 = seq_bin1;
        ss1 = seq_size2 ;
        ss2 = seq_size1 ;
    }

    long int it = 0;

    // Allocate memory and verify it has been allocated
    long int* xor = NULL;
    xor = calloc(ss1, sizeof(*xor));
    if (!xor)
        return printf("ERROR: xor_binary_array: cannot allocate memory.\n"), NULL;

    // Xor the two sequences values since its value is the same length.
    for (it = 0; it < ss2 - 1; it++)
        xor[it] = s1[it] ^ s2[it];

    // If one sequence is larger than the other, shift the last value of the smaller sequence toxor it with the other value.
    xor[it] = s1[it] ^ ((s2[it] << ((ss1 - ss2))));
    it++;

    // Values from the largest binary array are assigned to the xor result. (x^0 = x)
    for (it = ss2; it < ss1; it++)
        xor[it] = s1[it];

    return xor;
}

/**
 * Popcount a binary array sequence.
 * 
 * in : seq_bin : sequence in binary array format
 * in : seq_size : seq_bin length (number total of used bits)
 * out : bin_popcount : popcount of the seq : number of '1'
 * 
 * Iterates on seq_bin and for each value, adds its popcount to bin_popcount.
 */
int popcount_binary_array(const _Bool *seq_bin, const long int seq_size){
    int bin_popcount = 0;

    // Find the size of the binary array
    long int array_size = seq_size / int_SIZE;
    if(seq_size % int_SIZE != 0) array_size++;

    // Parse the binary array
    for (long int i = 0; i < array_size; ++i)
        bin_popcount += __builtin_popcount(seq_bin[i]);

    return bin_popcount;
}




/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/

//////////////// Convert to binary
/**
 * Convert a DNA base sequence to its binary array format.
 * 
 * in : dna_seq : DNA sequence in char array format
 * in : size : dna_seq length = number of nucleotides (= number of letters)
 * out : seq : DNA sequence in binary array format
 * 
 * Calls set_binary_array.
 */
long int* convert_to_binary(const char* dna_seq, const unsigned size){


    // Allocate memory and verify it has been allocated
    _Bool* seq_bin = NULL;
    seq_bin = calloc(size*2, sizeof(_Bool));
    if(!seq_bin)
        return printf("ERROR: set_binary_array: cannot allocate memory.\n"), NULL;

    int pos = 0;
    // Parse the DNA sequence, per nucleotides
    for (long int i = 0; i < size; ++i){
        // Set seq_bin bit values according to the nucleotide read
        switch(dna_seq[i]){
        case 'A': // A = 00
            break;
        case 'T': // T = 11
            seq_bin[i] = 1;
            seq_bin[i+1] = 1;
            break;
        case 'G': // G = 01
            seq_bin[i+1] = 1;
            break;
        case 'C': // C = 10
            seq_bin[i] = 1;
            break;
        case 'N': // N = A = 00
            break;
        case 'R': // R = A = 00
            break;
        case 'Y': // Y = C = 10
            seq_bin[i] = 1;
            break;
        case 'K': // K = G = 01
            seq_bin[i+1] = 1;
            break;
        case 'M': // M = A = 00
            break;
        case 'S': // S = C = 10
            seq_bin[i] = 1;
            break;
        case 'W': // W = A = 00
            break;
        case 'B': // B = C = 10
            seq_bin[i] = 1;
            break;
        case 'D': // D = A = 00
            break;
        case 'H': // H = G = 01
            seq_bin[i+1] = 1;
            break;
        case 'V': // V = A = 00
            break;
        }
        pos += 2;
    }
    return seq_bin;
}

//////////////// Convert binary aa to codon
/**
 * Convert a DNA sequence in binary array format to its DNA bases.
 * 
 * in : bin_dna_seq : DNA sequencDNA sequence in binary array format
 * in : size : total number total of used bits in the sequence bin_dna_seq
 * out : dna_seq : DNA sequence in char array format
 * 
 * For each pair of bits in bin_dna_seq, append to dna_seq its corresponding nucleotide.
 */
char* binary_to_dna(_Bool* bin_dna_seq, const unsigned size){
    if (size % 2 != 0) {
        printf("Error: binary_to_aa : wrong binary size (%d). Must be odd.\nExit.\n",size);
        return NULL;
    }

    //Allocate memory and verify it has been allocated
    char* dna_seq = calloc((size / 2) + 1, sizeof(*dna_seq));
    if(!dna_seq)
        return printf("ERROR: binary_to_dna: cannot allocate memory.\n"), NULL;

    int j = 0;
    //Parse the binary array, two bits per iteration
    for (unsigned i = 0; i < size; i += 2){
        // nucleotides = A, T, G, C
        int nucl1 = bin_dna_seq[i];
        int nucl2 = bin_dna_seq[i + 1];

        if (nucl1 == 0 && nucl2 == 0)
            // 00 = A/N
            dna_seq[j] = 'A';
        else if (nucl1 == 1 && nucl2 == 1)
            // 11 = T
            dna_seq[j] = 'T';
        else if (nucl1 == 1 && nucl2 == 0)
            // 10 = C
            dna_seq[j] = 'C';
        else if (nucl1 == 0 && nucl2 == 1)
            // 01 = G
            dna_seq[j] = 'G';
        j++;
    }
    return dna_seq;
}

//////////////// Generating mRNA
/**
 * Convert a DNA sequence in binary array format to its mRNA sequence.
 * 
 * in : gene_seq : DNA sequence in binary array format
 * in : seq_size : number total of used bits in the sequence gene_seq
 * out : rna_seq : resulting mRNA sequence in char array format
 * Convert a binary DNA sequence to a string mRNA sequence
 * 
 * For each pair of bits in bin_dna_seq, append to dna_seq its corresponding nucleotide in mRNA. (T -> U)
 */
char* generating_mRNA(const _Bool* gene_seq, const long int seq_size) {
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Allocate memory and verify it has been allocated
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_size / 2) + 2);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;

    long stop = seq_size;
    // Parse the binary DNA sequence, two bits per iteration
    for (long int i = 0; i < stop; i += 2) {

        // nucleotides = A, U, G, C
        int nucl1 = gene_seq[i];
        int nucl2 = gene_seq[i + 1];

        if (nucl1 == 0 && nucl2 == 0)
            // 00 = A
            rna_seq[j] = 'A';
        else if (nucl1 == 1 && nucl2 == 1)
            // 11 = U
            rna_seq[j] = 'U';
        else if (nucl1 == 1 && nucl2 == 0)
            // 10 = C
            rna_seq[j] = 'C';
        else if (nucl1 == 0 && nucl2 == 1)
            // 01 = G
            rna_seq[j] = 'G';
        else
            return printf("ERROR: generating_mRNA: invalid value in DNA sequence\n"), NULL;
        j++;
    }
    rna_seq[j] = '\0';
    return rna_seq;
}

//////////////// Detecting genes 
/**
 * Detects genes in the mRNA sequence in binary array format and maps them.
 * 
 * in : gene : mRNA sequence in binary array format
 * in : gene_size : number total of used bits in the sequence gene
 * in : gene_map : gene mapping struct
 * out : void
 * 
 * Iterates in the gene, and for each packet of 6 binary bits (corresponding to a nucleotide), searches for a start codon (AUG).
 * If a start codon is found, iterate until a stop codon is found (UAA, UAG or UGA).
 * If a stop codon is found, append to gene_map the gene length (from start to stop codon), its start position and stop one.
 * 
 * NB : The gene in binary array form can correspond to an mRNA or DNA sequence, since it is stored in the same way.
 */
void detecting_genes(const _Bool*gene, const long int gene_size, gene_map_t* gene_map) {
    gene_map->genes_counter = 0;

    // Check if memory ever have been allocated and allocate it if not
    if(!gene_map->gene_start || !gene_map->gene_end){
        gene_map->gene_start = malloc(sizeof(*gene_map->gene_start) * MAX_GENES);
        gene_map->gene_end = malloc(sizeof(*gene_map->gene_end) * MAX_GENES);

        if (!gene_map->gene_start || !gene_map->gene_end){
            printf("ERROR: detecting_genes: cannot allocate memory\n");
            return;
        }
    }

    int start_pos = -1;

    long int i = 0;

    //Parse the binary array, and find all the start and stop codons
    while ((i + 6) <= gene_size) {
        // Each nucleotides can be A, U, G or C
        int nucl1 = gene[i];
        int nucl2 = gene[i + 1];
        int nucl3 = gene[i + 2];
        int nucl4 = gene[i + 3];
        int nucl5 = gene[i + 4];
        int nucl6 = gene[i + 5];

        //If a start pos and a stop pos doesn't exist, search for AUG
        if (nucl1 == 0 && nucl2 == 0 && nucl3 == 1 
            && nucl4 == 1 && nucl5 == 0 && nucl6 == 1) {
        //if AUG, it's the start of a gene
            start_pos = i;
            i += 6;
        }
        else{

            if (start_pos != -1 ) {
                //if a start pos exists , search for UAA / UAG / UGA
                if ((nucl1 == 1 && nucl2 == 1 && nucl3 == 0) 
                    && ((nucl4 == 0 && nucl5 == 0 && nucl6 == 0)
                        || (nucl4 == 0 && nucl5 == 0 && nucl6 == 1)
                        || (nucl4 == 1 && nucl5 == 0 && nucl6 == 0))) {
                   //It's the end of a gene          
                   //If a start pos and an stop pos has been found, a gene exists so we save it in the struc
                    gene_map->gene_start[gene_map->genes_counter] = start_pos;
                    gene_map->gene_end[gene_map->genes_counter] = i+5;

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

/**
 * Retrives amino acid chains in a mRNA sequence in binary array format.
 * 
 * in : gene_seq : DNA sequence in binary array format
 * in : seq_size : gene_seq length (number total of used bits)
 * out : aa_seq : char array of proteins symbols.
 * 
 * The program parses the mRNA sequence, verify its length (seq_size).
 * Then iterates on gene_seq and for each packet of 6 binary bits (corresponding to a nucleotide), append to aa_seq its corresponding protein symbol.
 * 
 * NB : The gene in binary array form can correspond to an mRNA or DNA sequence, since it is stored in the same way.
*/
char* generating_amino_acid_chain(const _Bool*gene_seq, const long int seq_size) {
    long int codon_size = 6;
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;
    if(seq_size % 3 != 0)
        return NULL;

    // Allocate memory and verify it has been allocated
    char* aa_seq = NULL;
     aa_seq = malloc(sizeof(*aa_seq) * (seq_size / codon_size) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned temp = 0;

    long size = seq_size;

    //Parse the binary array, six bits by six (to parse three nucleotides per three)
    for (long int i = 0; i < size; i += codon_size) {
        // The hash functions, takes the 6 bits, and transform the array into an integer.
        // The integer first char is a 2, for hash generation purposes.
        int hash = 2;
        for(long int k = i; k < i + codon_size; k++){
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
            return NULL;
        }
    
        temp++;
    }

    aa_seq[temp] = '\0';
    return aa_seq;
}


/**
 * Detects probable mutation areas.
 * 
 * in : gene_seq : DNA sequence in binary array format
 * in : size_sequence : gene_seq length (number total of used bits)
 * in : mut_m : map of the possible mutation's areas
 * out : void
 * 
 * The algorithm scans a gene sequence and locates the high frequency of GC DNA bases in the sequence.
 * Must be at least 1/5th of the gene length
 * Precondition: gene_seq is of size size_sequence.
 * 
 * NB : The gene in binary array form can correspond to an mRNA or DNA sequence, since it is stored in the same way.
 */
void detecting_mutations(const _Bool*gene_seq, const long int size_sequence,
                         mutation_map mut_m) {
    long int detect_mut = 0;  //Counting size of GC sequence
    unsigned short tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones

    long size = size_sequence;
    //Parse the binary array, from the 'start_pos' bit to the end
    for (long int i = 0; i < size; i += 2) {

        // each nucleotides can be  A, U, G or C
        int nucl1 = gene_seq[i];
        int nucl2 = gene_seq[i + 1];

        //Increment detect_mut if find a C or G nucl
        if (((nucl1 == 0) && (nucl2 == 1)) ||
            ((nucl1 == 1) && (nucl2 == 0))) {
            if(detect_mut == 0){tmp_start_mut = i;}
            detect_mut+=2;
        }
        //Put detect_mut to 0 if find a A or T nucl
        else {
            //Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= (size_sequence / 5)) {
                mut_m.start_mut[cmp] = tmp_start_mut;
                mut_m.end_mut[cmp] = (i);
                mut_m.size[cmp] = detect_mut-1;
                cmp++;
            }
            detect_mut = 0;
        }
    }
    //Check if ending sequence is a probable mutation zone
    if (detect_mut >= (size_sequence / 5)) {
        mut_m.start_mut[cmp] = tmp_start_mut;
        mut_m.end_mut[cmp] = size_sequence;
        mut_m.size[cmp] = detect_mut-1;
    }
}

/**
 * Calculates the matching score of two binary array sequences.
 * 
 * in : seq1 : first sequence in binary
 * in : sequence_size1 : number total of used bits in the sequence seq1
 * in : seq2 : second sequence in binary
 * in : sequence_size2 : number total of used bits in the sequence seq2
 * out : float : 
 * 
 * The algorithms runs the hamming distance between two binary sequences, and return their matching score percentage
*/
float calculating_matching_score(const _Bool*seq1,const int seq_size1,
                                 const _Bool*seq2,const int seq_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    // First step: apply the xor operation between both arrays 

    long int *xor = NULL;
    xor = xor_binary_array(seq1, seq_size1, seq2, seq_size2);
    // xor_size = max size between 'seq_size1' and 'seq_size2'
    int xor_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;

    //Second step: count the number of bits whose value is 1 on the result
    int pop = popcount_binary_array(xor, xor_size);

    //Last step: compute the percentage
    float y = ((float)pop * 100.0) / (float)xor_size;
    return 100.0 - y;
}