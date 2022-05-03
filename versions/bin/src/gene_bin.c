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
    int new_pos = pos / 64;
    return __rolq(seq_bin[new_pos], pos+1) & 1;
}


/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/

//////////////// Convert to binary
/**
 * Convert a DNA base sequence to its binary array format.
 * 
 * in : seq_char : DNA sequence in char array format
 * in : seq_char_size : dna_seq length = number of nucleotides (= number of letters)
 * in/out : seq_bin : DNA sequence in binary array format
 * 
 * Set each int64 element of seq_bin from bit values according to the nucleotide read.
 * The non-ACGT nucleotides corresponding to several possible nucleotides are arbitrarily defined.
 */
void convert_to_binary(mm_array_t *seq_bin, const uint64_t seq_bin_size, const char* seq_char, const uint64_t seq_char_size)
{
    uint64_t k = 0;

    for(uint64_t j = 0; j < seq_bin_size - 1; ++j)
    {
        uint64_t p1 = 0, p2 = 0, p3 = 0, p4 = 0;
#ifdef __AVX512__
        uint64_t p5 = 0, p6 = 0, p7 = 0, p8 = 0;
#endif
        int i = 0;
        for(i = i; i < 32; ++i)
        {   
            p1 <<= 1;
            p1 = ((p1 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 64; ++i)
        {   
            p2 <<= 1;
            p2 = ((p2 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 96; ++i)
        {   
            p3 <<= 1;
            p3 = ((p3 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 128; ++i)
        {   
            p4 <<= 1;
            p4 = ((p4 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

#ifdef __AVX512__
        for(i = i; i < 160; ++i)
        {   
            p5 <<= 1;
            p5 = ((p5 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 192; ++i)
        {   
            p6 <<= 1;
            p6 = ((p6 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 224; ++i)
        {   
            p7 <<= 1;
            p7 = ((p7 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        for(i = i; i < 256; ++i)
        {   
            p8 <<= 1;
            p8 = ((p8 + L[seq_char[i+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[i+k]-OFFSET_TABLE][1];
        }

        seq_bin[j].reg = _mm512_setr_epi64x(p8, p7, p6, p5, p4, p3, p2, p1);
        k += 256;
#else
        seq_bin[j].reg = _mm256_setr_epi64x(p4, p3, p2, p1);
        k += 128;
#endif
    }

    uint64_t rest = seq_char_size - k;
#ifdef __AVX512__
    uint64_t p[8] = {0,0,0,0,0,0,0,0};
#else
    uint64_t p[4] = {0,0,0,0};
#endif
    int i;
    int r = rest >> 5;
    int m = rest & 31;
    for(i = 0; i < r; ++i)
    {   
        for(int j = 0; j < 32; ++j)
        {
            p[i] <<= 1;
            p[i] = ((p[i] + L[seq_char[j+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[j+k]-OFFSET_TABLE][1];
        }
        k += 32;
    }

    for(int j = 0; j < m; ++j)
    {
        p[i] <<= 1;
        p[i] = ((p[i] + L[seq_char[j+k]-OFFSET_TABLE][0]) << 1) + L[seq_char[j+k]-OFFSET_TABLE][1];
    }
    p[i] <<= 64 - 2*m;

#ifdef __AVX512__
    seq_bin[seq_bin_size - 1].reg = _mm512_setr_epi64x(p[7], p[6], p[5], p[4], p[3], p[2], p[1], p[0]);
#else
    seq_bin[seq_bin_size - 1].reg = _mm256_setr_epi64x(p[3], p[2], p[1], p[0]);
#endif
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

    int j = 0;
    //Parse the binary array, two bits per iteration
    for (unsigned i = 0; i < size; i += 2){
        // nucleotides = A, T, G, C
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
 * Convert a DNA sequence in binary array format to its mRNA sequence.
 * 
 * in : gene_seq : DNA sequence in binary array format
 * in : seq_size : number total of used bits in the sequence gene_seq
 * out : rna_seq : resulting mRNA sequence in char array format
 * Convert a binary DNA sequence to a string mRNA sequence
 * 
 * For each pair of bits in bin_dna_seq, append to dna_seq its corresponding nucleotide in mRNA. (T -> U)
 */

//TODO: faire à partir du string et non depuis à partir du binaire
char* generating_mRNA(const long int* gene_seq, const long start_pos, const long int seq_size) {
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Allocate memory and verify it has been allocated
    char* rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * (seq_size / 2) + 2);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    int j = 0;

    long stop = seq_size+start_pos;
    // Parse the binary DNA sequence, two bits per iteration
    for (long int i = start_pos; i < stop; i += 2) {

        // nucleotides = A, U, G, C
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
    rna_seq[j] = '\0';
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
void detecting_genes(const long int *gene, const long int gene_size, gene_map_t* gene_map)
{
    gene_map->genes_counter = 0;

    if (!gene_map->gene_start || !gene_map->gene_end)
    {
        printf("ERROR: detecting_genes: cannot allocate memory\n");
        return;
    }

    int start_pos = -1;

    long int i = 0;

    //Parse the binary array, and find all the start and stop codons
    while ((i + 6) <= gene_size)
    {
        int nucl = 0;
        // Each nucleotides can be A, U, G or C
        nucl = (nucl + get_binary_value(gene, i)) << 1;
        nucl = (nucl + get_binary_value(gene, i + 1)) << 1;
        nucl = (nucl + get_binary_value(gene, i + 2)) << 1;
        nucl = (nucl + get_binary_value(gene, i + 3)) << 1;
        nucl = (nucl + get_binary_value(gene, i + 4)) << 1;
        nucl = (nucl + get_binary_value(gene, i + 5));

        //If a start pos and a stop pos doesn't exist, search for AUG
        if (nucl == 13)
        {
        //if AUG, it's the start of a gene
            start_pos = i;
            i += 6;
        }
        else
        {
            //if a start pos exists , search for UAA / UAG / UGA
            if (start_pos != -1 && (nucl == 48 || nucl == 49 || nucl == 52))
            {
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
    // if(seq_size % 3 != 0)
    //     return NULL;

    // Allocate memory and verify it has been allocated
    char* aa_seq = malloc(sizeof(char) * (seq_size / codon_size) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned temp = 0;

    long size = start_pos+seq_size;

    //Parse the binary array, six bits by six (to parse three nucleotides per three)
    for (long int i = start_pos; i < size; i += codon_size) {
        // The hash functions, takes the 6 bits, and transform the array into an integer.
        // The integer first char is a 2, for hash generation purposes.
        int hash = 2;
        for(long int k = i; k < i + codon_size; k++){
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
            return NULL;
        }
    
        temp++;
    }

    aa_seq[temp] = '\0';
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
    //Parse the binary array, from the 'start_pos' bit to the end
    for (long int i = start_pos; i < size; i += 2) {

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

    int size = (seq_size1 > seq_size2) ? seq_size1 : seq_size2;

    // nombre d'éléments dans tableau
    int nb = size/ int_SIZE;
    if(size % int_SIZE != 0)
        nb++;

    long int x = 0, pop = 0;
    for(int i = 0; i < nb; ++i)
    {
        x = seq1[i] ^ seq2[i];
        pop += __builtin_popcountl(x);
    }

    return 100.0 - (float)(pop * 100.0 / (size));
}
