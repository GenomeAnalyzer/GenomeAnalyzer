#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "gene_bin.h"


/***************************************/
/********** BINARIES FUNCTION **********/
/***************************************/


/**
 * in : seq_bin : "binary array"
 * in : pos : position of the wanted bit (0 <= pos < size_seq * 2)
 * out : int : the wanted bit
 */
int get_binary_value(const unsigned int *seq_bin, const int pos){
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;
    return seq_bin[new_pos] & (1 << (pos - new_pos)) ? 1 : 0;
}

/**
 * in : seq_bin : "binary array"
 * in : pos : position of the bit to be replaced (0 <= pos < size_seq * 2)
 * in : value : the bit to replace
 * out : seq_bin : "binary array" modified
 */
unsigned int* change_binary_value(unsigned int *seq_bin, const int pos, const int value){
    // value = 1 or 0
    // ~(1000) = (0001)
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;

    if (value)
        seq_bin[new_pos] |= (1 << (pos - new_pos));
    else
        seq_bin[new_pos] &= ~(1 << (pos - new_pos));
    return seq_bin;
}

/**
 * in : seq_char : the DNA seq in it's char* mode
 * in : seq_size : size of the array 'seq_char' = number of nucleotides
 * out : seq_bin : "binary array" filled
 */
unsigned int* set_binary_array(const char *seq_char, const unsigned seq_size){
    int seq_bin_size = 2 * seq_size;
    int nb = seq_bin_size / int_SIZE;

    if(seq_bin_size % int_SIZE != 0)
        nb++;
    unsigned int* seq_bin = NULL;
    seq_bin = malloc(sizeof(*seq_bin) * nb);

    for (int i = 0; i < seq_bin_size * 2; ++i)
        change_binary_value(seq_bin, i, 0);

    int pos = 0;
    for (int i = 0; i < seq_size; ++i)
    {
        switch(seq_char[i]){
        case 'A':
            // A = 00
            change_binary_value(seq_bin, pos, 0);
            change_binary_value(seq_bin, pos + 1, 0);
            break;
        case 'T':
            // T = 11
            change_binary_value(seq_bin, pos, 1);
            change_binary_value(seq_bin, pos + 1, 1);
            break;
        case 'G':
            // G = 01
            change_binary_value(seq_bin, pos, 0);
            change_binary_value(seq_bin, pos + 1, 1);
            break;
        case 'C':
            // C = 10
            change_binary_value(seq_bin, pos, 1);
            change_binary_value(seq_bin, pos + 1, 0);
            break;
        }
        pos += 2;
    }
    return seq_bin;
}

/**
 * in : seq_bin1 : first sequence to xor: "binary array"
 * in : seq_size1 : number total of used bits in the sequence seq_bin1
 * in : seq_bin2 : second sequence to xor: "binary array"
 * in : seq_size2 : number total of used bits in the sequence seq_bin2
 * out : xor : "binary array" which is the result of the xor operation beteween seq1 and seq2
 */
unsigned int* xor_binary_array(const unsigned int *seq_bin1, const unsigned seq_size1,
                                const unsigned int *seq_bin2, const unsigned seq_size2){

    // "nb" = number of 'int' in the arrays 'seq1' and 'seq2' 
    //      = size of the 'unsigned int*' arrays
    int nb_seq1 = seq_size1 / int_SIZE;
    if(seq_size1 % int_SIZE != 0)   nb_seq1++;
    int nb_seq2 = seq_size2 / int_SIZE;
    if(seq_size2 % int_SIZE != 0)   nb_seq2++;

    int max_nb = nb_seq1 >= nb_seq2 ? nb_seq1 : nb_seq2;

    int max_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;
    int min_size = seq_size1 >= seq_size2 ? seq_size2 : seq_size1;
    int diff_size = max_size - min_size;

    unsigned int *xor = NULL;
    xor = malloc(sizeof(*xor) * max_nb);


    // If the sequences don't have the same size, do (with x = 1 or 0):

    //  xxxxxxxxx
    // ^
    //  000xxxxxx
    // -----------
    //  xxxxxxxxx

    // And, 0 ^ x = x

    int bit = 0;
    if(max_size == seq_size1){
        for (int i = 0; i < diff_size; ++i){
           bit = get_binary_value(seq_bin1, i);
           xor = change_binary_value(xor, i, bit);
        }
        for (int i = diff_size; i < max_size; ++i){
           bit = get_binary_value(seq_bin1, i) ^ get_binary_value(seq_bin2, i - diff_size);
           xor = change_binary_value(xor, i, bit);
        }
    }
    else{
        for (int i = 0; i < diff_size; ++i){
           bit = get_binary_value(seq_bin2, i);
           xor = change_binary_value(xor, i, bit);
        }
        for (int i = diff_size; i < max_size; ++i){
           bit = get_binary_value(seq_bin1, i - diff_size) ^ get_binary_value(seq_bin2, i);
           xor = change_binary_value(xor, i, bit);
        }
    }
    return xor;
}

/**
 * in : seq_bin : "binary array"
 * in : seq_size : number total of used bits in the sequence
 * out : s : popcount of the seq : number of '1'
 */
int popcount_binary_array(const unsigned int *seq_bin, const unsigned int seq_size){
    int s = 0;

    int nb = seq_size / int_SIZE;
    if(seq_size % int_SIZE != 0)    nb++;

    for (int i = 0; i < nb; ++i)
        s += __builtin_popcount(seq_bin[i]);

    return s;
}

/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/

//////////////// Convert to binary
/**
 * in : dna_seq : array of char
 * in : size : size of dna_seq = number of nucleotides (= number of letter)
 * out : seq : array of int
 * Convert a char DNA sequence to its binary sequence
 */
unsigned int* convert_to_binary(const char* dna_seq, const unsigned size){
    return set_binary_array(dna_seq, size);
}

//////////////// Convert binary aa to codon
/**
 * in : bin_dna_seq : unsigned int array
 * in : size : number total of used bits in bin_dna_seq
 * out : dna_seq : array of binary
 * Convert a binary sequence to its DNA bases
 */
char* binary_to_dna(unsigned int* bin_dna_seq, const unsigned size){
    if (size % 2 != 0) {
        printf("Error: binary_to_aa : wrong binary size (%d). Must be odd.\nExit.\n",size);
        return NULL;
    }

    char* dna_seq = malloc(sizeof(*dna_seq) * (size / 2) + 1);

    int j = 0;
    for (unsigned i = 0; i < size; i += 2){
        // nucleotides = A, T/U, G, C
        int nucl1 = get_binary_value(bin_dna_seq, i);
        int nucl2 = get_binary_value(bin_dna_seq, i + 1);

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
 * in : gene_seq : unsigned int array : "binary array"
 * in : seq_size : size of gene_seq : number of used bits in gene_seq
 * out : rna_seq : array of char
 * Convert a binary DNA sequence to a string mRNA sequence
 */
char* generating_mRNA(const unsigned int* gene_seq, const unsigned int seq_size) {
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

        // nucleotides = A, T/U, G, C
        int nucl1 = get_binary_value(gene_seq, i);
        int nucl2 = get_binary_value(gene_seq, i + 1);

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
    return rna_seq;
}

//////////////// Detecting genes 
/**
* in : gene : unsigned int array : "binary array"
* in : gene_size : size of gene : number of used bits in gene
* in : gene_map : struct to map the genes
* out : void
* Detect if a gene exists in the sequence and insert it in the structure
*/
void detecting_genes(const unsigned int *gene, const unsigned int gene_size, gene_map_t* gene_map) {
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
        // nucleotides = A, T/U, G, C
        int nucl1 = get_binary_value(gene, i);
        int nucl2 = get_binary_value(gene, i + 1);
        int nucl3 = get_binary_value(gene, i + 2);
        int nucl4 = get_binary_value(gene, i + 3);
        int nucl5 = get_binary_value(gene, i + 4);
        int nucl6 = get_binary_value(gene, i + 5);

        //If a start pos and a stop pos doesn't exist, search for AUG
        // if (!(gene[i%32] & ( 1 << (i%32) ))
        if (nucl1 == 0 && nucl2 == 0 && nucl3 == 1 
            && nucl4 == 1 && nucl5 == 0 && nucl6 == 1) {
        //if atc, it's the start of a gene
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

/**
 * in : gene_seq : unsigned int array : "binary array"
 * in : seq_size : size of gene_seq : number of used bits in gene_seq
 * out : char* : protein in symbols
 * The program parses the mRNA sequence, verify its length and if the first codon is a START codon.
*/
char* generating_amino_acid_chain(const unsigned int *gene_seq, const unsigned int seq_size) {
    short codon_size = 6;
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;

    // Create and check the output
    char* aa_seq = NULL;
    aa_seq = malloc(sizeof(*aa_seq) * (seq_size / codon_size) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned temp = 0;

    for (unsigned int i = 0; i < seq_size; i += codon_size) {
        // The hash functions, takes the 6 bits, and transform the array into an integer.
        // The integer first char is a 2, for hash generation purposes.
        int hash = 2;
        for(int k = i; k < i + codon_size; k++){
            hash = 10 * hash + get_binary_value(gene_seq, k);
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
    return aa_seq;
}


/**
* in : gene_seq : sequence of the gene : "binary array"
* in : size_sequence : size of the sequence : number of used bits in gene_seq
* in : mut_m : map of the possible mutation's areas
* out : boolean
* The algorithm runs through a gene sequence and detects if there is a sequence with a high GC frequency
* (at least 1/5th of the gene sequence's size) it returns true, else it returns false.
* Precondition: gene_seq is of size size_sequence.
*/
void detecting_mutations(const unsigned int *gene_seq, const unsigned int size_sequence,
                         mutation_map mut_m) {
    unsigned long detect_mut = 0;  //Counting size of GC sequence
    unsigned short tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones

    //Read the sequence
    for (unsigned long i = 0; i < size_sequence; i += 2) {

        // nucleotides = A, T/U, G, C
        int nucl1 = get_binary_value(gene_seq, i);
        int nucl2 = get_binary_value(gene_seq, i + 1);

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

/**
 * in : seq1 : first sequence in binary
 * in : sequence_size1 : number total of used bits in the sequence seq1
 * in : seq2 : second sequence in binary
 * in : sequence_size2 : number total of used bits in the sequence seq2
 * out : float : matching score in %
 * The algorithms runs the hamming distance between two binary sequences, and return their matching score in %
*/
float calculating_matching_score(const unsigned int *seq1, const int seq_size1,
                                 const unsigned int *seq2, const int seq_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    unsigned int *xor = NULL;
    xor = xor_binary_array(seq1, seq_size1, seq2, seq_size2);

    // xor_size = max size
    int xor_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;

    int pop = popcount_binary_array(xor, xor_size);

    float y = ((float)pop * 100.0) / (float)xor_size;
    return 100.0 - y;
}
