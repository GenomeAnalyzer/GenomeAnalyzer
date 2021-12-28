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
 * in : seq_size : size of the array 'seq_char' = number of nitrogenous bases
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

/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/