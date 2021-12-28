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



/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/