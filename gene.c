#include "gene.h"
#include <stdio.h>
#include <stdlib.h>

//A   : 00
//T/U : 11
//C   : 10
//G   : 01

// UAA 110000
// UAG 110001
// UGA 110100

// ATC 001110

//////////////// Generating mRNA
void generating_mRNA(){
    // for (unsigned long long i = 0; i < seq_len; i++)
    //     if (seq[i] == 'T')
    //         seq[i] = 'U';
}


//////////////// Algo naif
void detecting_genes(unsigned int gene[],gene_map_t *gene_map){


    //struct gene_map_s gene_map;
    gene_map->genes_counter = 0;

    



    int start_pos = -1;
    int stop_pos = -1;

    int i = 0;

    while((i+6) < 1000)
    {
        

    if (start_pos == -1 && stop_pos == -1)
    {

        // if (!(gene[i%32] & ( 1 << (i%32) ))
        if(gene[i] == 0 && gene[i+1] == 0 && gene[i+2] == 1 && gene[i+3] == 1 && gene[i+4] == 1 && gene[i+5] == 0)
        {
            start_pos = i;
        }
    }else
    if(start_pos != -1 && stop_pos == -1)
    {
         if(((gene[i] == 1 && gene[i+1] == 1 && gene[i+2] == 0) && (gene[i+3] == 0 && gene[i+4] == 0 && gene[i+5] == 0) )||(gene[i+3] == 0 && gene[i+4] == 0 && gene[i+5] == 1)||(gene[i+3] == 1 && gene[i+4] == 0 && gene[i+5] == 0))
         {
            stop_pos = i;  
         }     
    }else{

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
void generating_amino_acid_chain(){

}

//////////////// Detecting probable mutation zones
void detecting_mutations(){

}

//////////////// Calculating the matching score of two sequences
void calculating_matching_score(){

}

//////////////// Hamming calculation
void hamming(){
    // for (i = 0; i < min(lengt(seq0), lenghth(seq1)); i++)
    //     d += popcount(seq0[i] xor seq1[i]);
}