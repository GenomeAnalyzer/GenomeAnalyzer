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


/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/