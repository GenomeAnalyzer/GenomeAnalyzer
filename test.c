#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.c"
#include <string.h>

int main(int argc, char* argv[]){
    // char* seq = "AAAAA";
    // short* res = compress(5, seq);
    // for(int i =0; i< 10; i++){
    //     printf("%d_",res[i]);
    // }

    // int tab[2] = {1, 2};
    // int taby[2] = {1, 2};
    // // printf("%p, %p\n", tab, taby);
    // // printf("%d, %d\n", sizeof(tab), sizeof(taby));
    // // printf("%d\n", memcmp(tab, taby, sizeof(tab)));
    // // printf("%ld_%ld\n", *tab, *taby);
    // if (!memcmp(tab, taby, sizeof(tab))) printf("yey !");
    // else printf("nono");

    // char gene [] = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT";
    // "001111000000010111111100110010101111101010000101110000100000001010000010100000101111111001001110111011110111000100111011"
    short seq [] = { 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1 };
    short seq_size = 120;

    gene_map_t gene_map;
    gene_map.genes_counter = 3;
    gene_map.gene_start = malloc(sizeof(short)*3);
    gene_map.gene_start[0] = 0;
    gene_map.gene_start[1] = 12;
    gene_map.gene_start[2] = 96;
    gene_map.gene_end = malloc(sizeof(short)*3);
    gene_map.gene_end[0] = 12;
    gene_map.gene_end[1] = 24;
    gene_map.gene_end[2] = 120;

    short codons_count = 2;
    codon codons[2];
    // test
    short tab[6] = { 0,0,1,1,1,1 };
    memcpy(codons[0].codon, tab, sizeof(short) * 6);
    short tab2[6] = { 0,0,0,0,0,0 };
    memcpy(codons[1].codon, tab2, sizeof(short) * 6);

    codons[0].symbol = 'I';
    codons[1].symbol = 'K';

    

    char** proteines = malloc(sizeof(gene_map.genes_counter));
    
    // char* proteines[seq_size/6];
    // printf("%p\n", proteines);

    generating_amino_acid_chain(proteines, seq_size, seq, gene_map, codons_count, codons);

    printf("%p\n", proteines);
    for (int i  = 0; i < gene_map.genes_counter; i++)
        printf("proteines[%d] : %s\n", i, proteines[i]);
}