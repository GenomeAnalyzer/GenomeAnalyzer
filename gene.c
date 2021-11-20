#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.h"

//A   : 00
//T/U : 11
//C   : 10
//G   : 01

// UAA 110000
// UAG 110001
// UGA 110100

// ATC 001110

//TODO: optimize + exit when strange char
int* convert_to_binary(char* dna_seq)
{

    unsigned i = 0;
    unsigned temp = 0;
    
   // printf("\n%s FIN%d\n",dna_seq,sizeof(dna_seq));
    //dna_seq = "TGGGCCCA\n";


    int *seq = malloc(sizeof(int)*MAX_GENES); 
    while(dna_seq != "\n")
   //for(int c = 0; c <= sizeof(dna_seq); c ++ )
    {
        switch(dna_seq[i])
        {
            case 'A':
                seq[temp] = 0;
                seq[temp+1] = 0; 
                temp += 2;

                break;
            
            case 'T':
                seq[temp] = 1;
                seq[temp+1] = 1;
                temp += 2;

                break;
            
            case 'C':
                seq[temp] = 1;
                seq[temp+1] = 0;
                temp += 2;

                break;
            case 'G':
                seq[temp] = 0;
                seq[temp+1] = 1;
                temp += 2;

                break;
            case 'N':
                break;
            default:
                printf("Error: wrong letter in the sequence(%c).\nExit.\n",dna_seq[i]);
                    seq[temp] = NULL;

                return seq;

               // exit(-1);
        }
         i++; 

    }

    seq[temp] = NULL;
    return seq;

}

//////////////// Generating mRNA
/**
 * generating_mRNA - Convert a binary DNA sequence to a string mRNA sequence
 * dna_seq      array of binary
 * seq_len      size of dna_seq
 *
 * rna_seq      array of char
 */
char* generating_mRNA(int* dna_seq, int seq_len) {
    // Check the input argument
    if (!dna_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Create and check the output
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_len / 2));
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;
    // Parse the binary DNA sequence two by two
    for (int i = 0; i < seq_len; i += 2) {
        switch (dna_seq[i]) {
        case 0:
            if (dna_seq[i + 1] == 0)
                rna_seq[j] = 'A';
            else if (dna_seq[i + 1] == 1)
                rna_seq[j] = 'G';
            break;
        case 1:
            if (dna_seq[i + 1] == 0)
                rna_seq[j] = 'C';
            else if (dna_seq[i + 1] == 1)
                rna_seq[j] = 'U';
            break;
        default:
            return printf("ERROR: generating_mRNA: invalid value in DNA sequence\n"), NULL;
        }
        j++;
    }
    return rna_seq;
}

//////////////// Detecting genes 
/*
*
*   Input : unsigned int gene[] -> sequence of genes , gene_map -> Struct to map the genes
*   Output : void
*   Function: detect if a gene exists in the sequence and insert it in the structure
*
*/
void detecting_genes(unsigned int gene [], gene_map_t* gene_map) {
    //struct gene_map_s gene_map;
    gene_map->genes_counter = 0;

    int start_pos = -1;
    int stop_pos = -1;

    int i = 0;

    while ((i + 6) < 1000) {


        if (start_pos == -1 && stop_pos == -1) {

            //If a start pos and a stop pos doesn't exist, search for ATC
            // if (!(gene[i%32] & ( 1 << (i%32) ))
            if (gene[i] == 0 && gene[i + 1] == 0 && gene[i + 2] == 1 && gene[i + 3] == 1 && gene[i + 4] == 1 && gene[i + 5] == 0) {
                //if atc, it's the start of a gene
                start_pos = i;
            }
        }
        else
            if (start_pos != -1 && stop_pos == -1) {
                //if a start pos exists , search for UAA / UAG / UGA
                if (((gene[i] == 1 && gene[i + 1] == 1 && gene[i + 2] == 0) && (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 0)) || (gene[i + 3] == 0 && gene[i + 4] == 0 && gene[i + 5] == 1) || (gene[i + 3] == 1 && gene[i + 4] == 0 && gene[i + 5] == 0)) {
                    //It's the end of a gene          
                   //If a start pos and an stop pos has been found, a gene exists so we save it in the struc
                    gene_map->gene_start[gene_map->genes_counter] = start_pos;
                    gene_map->gene_end[gene_map->genes_counter] = i;

                    gene_map->genes_counter++;

                    start_pos = -1;
                    stop_pos = -1;
                }
            }

        i++;
    }
}

//////////////// Generating an amino acid chain (protein) 
/*
 * in : seq : original mRNA sequence.
 * out : char* : protein in symbols
 * The program parses the mRNA sequence, verify its length and if the first codon is a START codon.
*/

char* generating_amino_acid_chain(char* seq, int seq_size, int codons_count, codon codons []) {
    gene_map_t gene;
    char* protein = "";
    gene.genes_counter = 0;
    gene.gene_start = "AUG";
    gene.gene_end = NULL;

    if ((seq_size % 3) != 0)
        return printf("ERROR: generating_amino_acid_chain: Invalid sequence size\n"), NULL;

    for (int i = 0; i < seq_size; i += 3) {
        char cod = seq[i] + seq[i + 1] + seq[i + 2];
        if (i == 0 && cod != "AUG")
            return printf("ERROR: generating_amino_acid_chain: Sequence does not start with a start codon (AUG)\n"), NULL;
        for (int j = 0; j < codons_count; j++) {
            if (cod == codons[j].codon) {
                protein += codons[j].symbol;
                gene.genes_counter += 1;
                if (codons[j].symbol == "O") {
                    gene.gene_end = codons[j].codon;
                    return protein;
                }
                break;
            }
        }
    }
    return printf("ERROR: generating_amino_acid_chain: sequence has no END codon\n"), NULL;
}

//////////////// Detecting probable mutation zones
/*
* in : gene_seq : sequence of the gene
* in : size_sequence : size of the sequence
* out : boolean
* The algorithm runs through a gene sequence and detects if there is a sequence with a high GC frequency
* (at least 1/5th of the gene sequence's size) it returns true, else it returns false.
* Precondition: gene_seq is of size size_sequence.
*/
bool detecting_mutations(const unsigned int gene_seq [], const unsigned long long size_sequence) {
    unsigned long long detect_mut = 0;  //Counting size of GC sequence

    //Read the sequence
    for (unsigned long long i = 0; i < size_sequence; i += 2) {
        //Increment detect_mut if find a C or G gene
        if (((__builtin_popcount(gene_seq[i]) == 0) && (__builtin_popcount(gene_seq[i + 1]) == 1)) ||
            ((__builtin_popcount(gene_seq[i]) == 1) && (__builtin_popcount(gene_seq[i + 1]) == 0))) {
            detect_mut++;
        }
        //Put detect_mut to 0 if find a A or T gene
        else {
            //Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= ((size_sequence / 2) / 5)) {
                return true;
            }
            detect_mut = 0;
        }
    }
    //Check if ending sequence is a probable mutation zone
    if (detect_mut >= ((size_sequence / 2) / 5)) {
        return true;
    }
    return false;
}


//////////////// Calculating the matching score of two sequences
/*
 * in : seq1 : First sequence in binary
 * in : seq2 : Second sequence in binary
 * out : float : matching score in %
 * The algorithms runs the hamming distance between two binary sequences, and return their matching score in %
*/

float calculating_matching_score(int sequence_size, int seq1 [], int seq2 []) {
    // Check outside of function that seq1 and seq2 are the same size

    int total_size_sequence = 0;
    int total_hamming_distance = 0;

    for (int i = 0; i < sequence_size; i++) {
        int size = binary_size_count(seq1[i]);
        // Make check outside of function ?
        if (size != binary_size_count(seq2[i])) {
            printf("ERROR: generating_mRNA: wrong size sequence\n");
            return -1.0;
        }
        total_size_sequence += size;
        total_hamming_distance += hamming(seq1[i], seq2[i]);
    }

    return 100 * total_hamming_distance / total_size_sequence;
}

//////////////// Counting binary size
/*
 * in : b : binary number
 * out : int : count of bits
 * The algorithm return the count of bits the binary integer is.
*/
int binary_size_count(int b) {
    int size = 0;
    while (b) {
        size++;
        b >>= 1;
    }
    return size;
}

//////////////// Hamming calculation
/*
 * in : seq1 : First sequence in binary
 * in : seq2 : Second sequence in binary
 * out : int : count of bits
 * The algorithms calculate the hamming distance between two binary sequences
*/
int hamming(int seq1, int seq2) {
    int x = seq1 ^ seq2;
    return __builtin_popcount(x);
}