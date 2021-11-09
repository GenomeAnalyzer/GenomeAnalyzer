#include "gene.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int n = 1000;

unsigned int* alea_gen(){
    unsigned int* gene = malloc(n* sizeof(int));

    srand(getpid());


    for (int i = 0; i < n ; i++ )
    {   
        gene[i] = rand()%2; 
        printf(" %d",gene[i]);

    }
    return gene;
}


unsigned int* test_gen(){
     unsigned int* gene= calloc(sizeof(int*),MAX_GENES);
    

    //Example of a gene so the gene map register it 
    gene[0] = 0;
    gene[1] = 0;
    gene[2] = 1;
    gene[3] = 1;
    gene[4] = 1;
    gene[5] = 0;
    gene[6] = 1;
    gene[7] = 1;
    gene[8] = 0;
    gene[9]= 0;
    gene[10] = 0;
    gene[11] = 0; 

    return gene;  
}

int main(){

    //unsigned int* gene= test_gen();
    
    unsigned int* gene= alea_gen();


   gene_map_t gene_map;

    gene_map.genes_counter = 0;

    gene_map.gene_end = calloc(sizeof(int*),MAX_GENES);

    gene_map.gene_start = calloc(sizeof(int*),MAX_GENES);

    detecting_genes(gene,&gene_map);

    //Test print 

    printf("\n%llu \n",gene_map.genes_counter);

    return 0;

}