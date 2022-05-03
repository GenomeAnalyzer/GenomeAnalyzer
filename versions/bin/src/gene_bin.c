#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

#include "../headers/gene_bin.h"

/***************************************/
/********** BINARIES FUNCTION **********/
/***************************************/

/**
 * Retrieve one bit from the binary array sequence.
 * 
 * in : seq_bin : sequence in binary array format
 * in : pos : position of the requested bit (0 <= pos < size_seq * 2)
 * out : int : the requested bit
 * 
 * Shifts the seq_bin by pos.
 */
int get_binary_value(const long int *seq_bin, const int pos){
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;
    return seq_bin[new_pos] & ((long)1 << (pos - new_pos)) ? 1 : 0;
}

/**
 * Change one bit in the binary array sequence.
 * 
 * in : seq_bin : sequence in binary array format
 * in : pos : position of the bit to be replaced (0 <= pos < size_seq * 2)
 * in : value : new bit value
 * out : seq_bin : sequence in binary array format with the one bit changed
 * 
 * Set the bit value of seq_bin at pos position to value
 */
long int* change_binary_value(long int *seq_bin, const int pos, const int value){
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;

    if (value)
        seq_bin[new_pos] |= ((long)1 << (pos - new_pos));
    else
        seq_bin[new_pos] &= ~((long)1 << (pos - new_pos));
    return seq_bin;
}

/**
 * Convert a char formated DNA sequence to its binary array format.
 * 
 * in : seq_char : the DNA seq in it's char* mode
 * in : seq_size : size of the array 'seq_char' (number of nucleotides)
 * out : seq_bin : sequence in binary array format
 * 
 * Iterates over seq_char and sets seq_bin bit values according to the nucleotide read.
 * The non-ACGT nucleotides corresponding to several possible nucleotides are arbitrarily defined.
 */
long int* set_binary_array(const char *seq_char, const unsigned seq_size){
    // Number of bits needed to transform seq_char into a binary array.
    int seq_bin_size = 2 * seq_size;

    // Binary array new size
    int nb = seq_bin_size / int_SIZE;
    if(seq_bin_size % int_SIZE != 0)
        nb++;

    // Allocate memory and verify it has been allocated
    long int* seq_bin = NULL;
    seq_bin = calloc(nb, sizeof(long int));
    if(!seq_bin)
        return printf("ERROR: set_binary_array: cannot allocate memory.\n"), NULL;

    int pos = 0;
    long int i;

    // Parse the DNA sequence, per nucleotides
#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,64)
    for (i = 0; i < seq_size; ++i){

        //Default char is put to A to handle the warning : wrong size input
        //Add '00' bits in this case
        int c = 0;
        if(!seq_char[i]){
            printf("WARNING: set_binary_array: size input is different than the char sequence.\n");
        }
        else{
            c = seq_char[i] - 65;
        }
        //get the 2-bits value of char read
        int bit1, bit2;
        bit1 = L[c][0];
        bit2 = L[c][1];

        pos = 2*i;
        // Set seq_bin bit values according to the nucleotide read
        change_binary_value(seq_bin, pos, bit1);
        change_binary_value(seq_bin, pos + 1, bit2);
    }
}
    return seq_bin;
}

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
long int* xor_binary_array(long int* const seq_bin1, const unsigned seq_size1,
                           long int * const seq_bin2, const unsigned seq_size2){

    // size of the binary array type used
    long int intsize = int_SIZE + 1;

    long int* s1, * s2;
    long int ss1, ss2;
    long int sbs1, sbs2;

    // Find the greater binary array, and rename them
    if (seq_size1 >= seq_size2) { // if s1 is greater or equal than s2
        s1 = seq_bin1;
        s2 = seq_bin2;
        ss1 = seq_size1 / intsize + ((seq_size1 / intsize) % intsize != 0);
        sbs1 = seq_size1;
        ss2 = seq_size2 / intsize + ((seq_size2 / intsize) % intsize != 0);
        sbs2 = seq_size2;
    }
    else { // else if s2 is greater or equal than s1
        s1 = seq_bin2;
        s2 = seq_bin1;
        ss1 = seq_size2 / intsize + ((seq_size2 / intsize) % intsize != 0);
        sbs1 = seq_size2;
        ss2 = seq_size1 / intsize + ((seq_size1 / intsize) % intsize != 0);
        sbs2 = seq_size1;
    }

    if (ss1 == 0) ss1 = 1;
    if (ss2 == 0) ss2 = 1;

    long int it = 0;

    // Allocate memory and verify it has been allocated
    long int* xor = NULL;
    xor = calloc(ss1, sizeof(*xor));
    if (!xor)
        return printf("ERROR: xor_binary_array: cannot allocate memory.\n"), NULL;

    // Xor the two sequences values since its value is the same length.
#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,64)
    for (it = 0; it < ss2 - 1; it++)
        xor[it] = s1[it] ^ s2[it];

#pragma omp master
    {
    // If one sequence is larger than the other, shift the last value of the smaller sequence toxor it with the other value.
    xor[ss2-1] = s1[ss2-1] ^ ((s2[ss2-1] << ((sbs1 - sbs2) % intsize)));
    }
    // Values from the largest binary array are assigned to the xor result. (x^0 = x)
#pragma omp for schedule(static,64)
    for (it = ss2; it < ss1; it++)
        xor[it] = s1[it];
}
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
int popcount_binary_array(const long int *seq_bin, const long int seq_size){
    int bin_popcount = 0;

    // Find the size of the binary array
    long int array_size = seq_size / int_SIZE;
    if(seq_size % int_SIZE != 0) array_size++;

    long int i;

    // Parse the binary array
#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,64) reduction(+:bin_popcount)
    for (i = 0; i < array_size; ++i)
        bin_popcount += __builtin_popcount(seq_bin[i]);
}
    return bin_popcount;
}

/**
 * Retrieve a piece of the binary array sequence.
 * 
 * in : seq_bin : sequence in binary array format
 * in : size : size of the binary requested
 * out : piece_seq_bin : the requested part of seq_bin, from pos_start and size.
 * 
 * Iterates on seq_bin from pos_start, size times and gets for each iteration its binary value.
 */
long int* get_piece_binary_array(const long int* seq_bin, const long int pos_start, const long int size){
    //Find the size of the output
    long int array_size = (size) / int_SIZE + (size % int_SIZE != 0);
    
    // Allocate memory and verify it has been allocated
    long int *piece_seq_bin = NULL;
    piece_seq_bin = calloc(array_size, sizeof(*seq_bin));
    if(!piece_seq_bin) 
        return printf("ERROR: get_piece_binary_array: cannot allocate memory.\n"), NULL;

    // stop position.
    long int stop_pos = pos_start + size;
    long i;
    //Parse the binary array,
    //from the bit at 'pos_start' position to 'pos_stop' position
#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,64)
    for(i = pos_start ; i < stop_pos; i++){
        long tmp = i - pos_start;
        change_binary_value(piece_seq_bin, tmp, get_binary_value(seq_bin,i));
    }
}
    return piece_seq_bin;
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
    return set_binary_array(dna_seq, size);
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
char* binary_to_dna(long int* bin_dna_seq, const unsigned size){
    if (size % 2 != 0) {
        printf("Error: binary_to_aa : wrong binary size (%d). Must be odd.\nExit.\n",size);
        return NULL;
    }

    //Allocate memory and verify it has been allocated
    char* dna_seq = calloc((size / 2) + 1, sizeof(*dna_seq));
    if(!dna_seq)
        return printf("ERROR: binary_to_dna: cannot allocate memory.\n"), NULL;

    //Parse the binary array, two bits per iteration
    unsigned i;
#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,32)
    for (i = 0; i < size; i += 2){
        // nucleotides = A, T, G, C
        int nucl1 = get_binary_value(bin_dna_seq, i);
        int nucl2 = get_binary_value(bin_dna_seq, i + 1);

        //get the ASCII value according to bits value
        char value = bitstocharDNA[nucl2 + 2*nucl1];

        int index = i/2;

        dna_seq[index] = value;
    }
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
char* generating_mRNA(const long int* gene_seq, const long start_pos, const long int seq_size) {
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Allocate memory and verify it has been allocated
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_size / 2) + 2);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    long stop = seq_size+start_pos;

    // Parse the binary DNA sequence, two bits per iteration
    long int i;

#pragma omp parallel default(shared)
{
#pragma omp for schedule(static,32)
    for (i = start_pos; i < stop; i += 2) {

        // nucleotides = A, U, G, C
        int nucl1 = get_binary_value(gene_seq, i);
        int nucl2 = get_binary_value(gene_seq, i + 1);

        long index = (i+1 - start_pos)/2;

        //get the ASCII value according to bits value
        char value = bitstocharmRNA[nucl2 + 2*nucl1];
        rna_seq[index] = value;
    }
}
    rna_seq[stop/2] = '\0';
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
void detecting_genes(const long int *gene, const long int gene_size, gene_map_t* gene_map) {
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
        int nucl1 = get_binary_value(gene, i);
        int nucl2 = get_binary_value(gene, i + 1);
        int nucl3 = get_binary_value(gene, i + 2);
        int nucl4 = get_binary_value(gene, i + 3);
        int nucl5 = get_binary_value(gene, i + 4);
        int nucl6 = get_binary_value(gene, i + 5);

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
char* generating_amino_acid_chain(const long int *gene_seq, const long int start_pos, const long int seq_size) {
    long int codon_size = 6;
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;
    if(seq_size % 3 != 0)
        return NULL;

    // Allocate memory and verify it has been allocated
    char* aa_seq = calloc(sizeof(char), sizeof(*aa_seq) * (seq_size / codon_size) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    long size = start_pos+seq_size;

    long int i,k;

    //Parse the binary array, six bits by six (to parse three nucleotides per three)
#pragma omp parallel default(shared) private(k)
{
#pragma omp for schedule(static,64/codon_size)
    for (i = start_pos; i < size; i += codon_size) {

        unsigned temp = (i - start_pos)/codon_size;

        //Get the decimal value of the 6 bits
        int tmp = 0;
        int pow_bit = 5;
        for(k = i; k < i + codon_size; k++){
            int get_bin = get_binary_value(gene_seq, k);
            tmp += get_bin << pow_bit;
            pow_bit--;
        }
        
        //Get the corresponding protein from the lookup table
        aa_seq[temp] = LUT[tmp];
    }
}
    aa_seq[size/codon_size] = '\0';
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
void detecting_mutations(const long int *gene_seq, const long int start_pos, const long int size_sequence,
                         mutation_map mut_m) {
    long int detect_mut = 0;  //Counting size of GC sequence
    unsigned short tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones

    long size = start_pos + size_sequence;
    long int i;

    //Parse the binary array, from the 'start_pos' bit to the end
    for (i = start_pos; i < size; i += 2) {

        // each nucleotides can be  A, U, G or C
        int nucl1 = get_binary_value(gene_seq, i);
        int nucl2 = get_binary_value(gene_seq, i + 1);

        //Increment detect_mut if find a C or G nucl
        if (((nucl1 == 0) && (nucl2 == 1)) ||
            ((nucl1 == 1) && (nucl2 == 0))) {
            if(detect_mut == 0){tmp_start_mut = i-start_pos;}
            detect_mut+=2;
        }
        //Put detect_mut to 0 if find a A or T nucl
        else {
            //Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= (size_sequence / 5)) {
                mut_m.start_mut[cmp] = tmp_start_mut;
                mut_m.end_mut[cmp] = (i)-start_pos;
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
float calculating_matching_score(const long int *seq1, long int start_pos1,const int seq_size1,
                                 const long int *seq2, long int start_pos2,const int seq_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    // First step: apply the xor operation between both arrays 
    long int *seq1tmp;
    seq1tmp = get_piece_binary_array(seq1, start_pos1, seq_size1);

    long int *seq2tmp;
    seq2tmp = get_piece_binary_array(seq2, start_pos2, seq_size2);    

    long int *xor = NULL;
    xor = xor_binary_array(seq1tmp, seq_size1, seq2tmp, seq_size2);
    // xor_size = max size between 'seq_size1' and 'seq_size2'
    int xor_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;

    //Second step: count the number of bits whose value is 1 on the result
    int pop = popcount_binary_array(xor, xor_size);

    //Last step: compute the percentage
    float y = ((float)pop * 100.0) / (float)xor_size;
    return 100.0 - y;
}
