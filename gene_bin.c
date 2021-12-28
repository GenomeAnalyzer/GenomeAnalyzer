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


/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/