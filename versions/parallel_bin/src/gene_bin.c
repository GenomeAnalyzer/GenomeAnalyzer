#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "../headers/gene_bin.h"
#include <sys/types.h>

#ifndef __USE_MISC
#define __USE_MISC
#endif
#include <dirent.h>

#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <nmmintrin.h>
#include <immintrin.h>

int output;

int rank;

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
int get_binary_value(const long int *seq_bin, const int pos)
{
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;
    return seq_bin[new_pos] & ((long)1 << (pos - new_pos)) ? 1 : 0;
}

/**
 * Change one bit in the binary array sequence.
 *
 * in : seq_bin : sequence in binary array format
 * in : pos : position of the bit to be replaced (0 <= pos < size_seq * 2)
 * in : value : new bit value
 * out : seq_bin : sequence in binary array format with the one bit changed
 *
 * Set the bit value of seq_bin at pos position to value
 */
long int *change_binary_value(long int *seq_bin, const int pos, const int value)
{
    int new_pos = pos > int_SIZE ? pos / int_SIZE : 0;

    if (value)
        seq_bin[new_pos] |= ((long)1 << (pos - new_pos));
    else
        seq_bin[new_pos] &= ~((long)1 << (pos - new_pos));
    return seq_bin;
}

/**
 * Convert a char formated DNA sequence to its binary array format.
 *
 * in : seq_char : the DNA seq in it's char* mode
 * in : seq_size : size of the array 'seq_char' (number of nucleotides)
 * out : seq_bin : sequence in binary array format
 *
 * Iterates over seq_char and sets seq_bin bit values according to the nucleotide read.
 * The non-ACGT nucleotides corresponding to several possible nucleotides are arbitrarily defined.
 */
long int *set_binary_array(const char *seq_char, const size_t seq_size)
{
    // Number of bits needed to transform seq_char into a binary array.
    int seq_bin_size = 2 * seq_size;

    // Binary array new size
    int nb = seq_bin_size / int_SIZE;
    if (seq_bin_size % int_SIZE != 0)
        nb++;

    // Allocate memory and verify it has been allocated
    long int *seq_bin = NULL;
    seq_bin = calloc(nb, sizeof(long int));
    if (!seq_bin)
        return printf("ERROR: set_binary_array: cannot allocate memory.\n"), NULL;

    int pos = 0;

    // int bit1, bit2, c;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static, 64)
        for (size_t i = 0; i < seq_size; ++i)
        {

            // Default char is put to A to handle the warning : wrong size input
            // Add '00' bits in this case
            int c = 0;
            if (!seq_char[i])
                printf("WARNING: set_binary_array: size input is different than the char sequence.\n");
            else
                c = seq_char[i] - 65;

            // get the 2-bits value of char read
            int bit1, bit2;
            bit1 = L[c][0];
            bit2 = L[c][1];

            pos = 2 * i;
            // Set seq_bin bit values according to the nucleotide read
            change_binary_value(seq_bin, pos, bit1);
            change_binary_value(seq_bin, pos + 1, bit2);
        }
    }
    return seq_bin;
}

/**
 * Xor two binary array sequences.
 *
 * in : seq_bin1 : first sequence in binary array format to xor
 * in : seq_size1 : seq_bin1 length (number total of used bits)
 * in : seq_bin2 : first sequence in binary array format to xor
 * in : seq_size2 : seq_bin2 length (number total of used bits)
 * out : xor : binary array sequence resulting from the xor operation between seq1 and seq2
 *
 * Iterates over sequences value per value.
 * Xor the two sequences values since its value is the same length.
 * If one sequence is larger than the other, shift the last value of the smaller sequence to xor it with the other value.
 * Values from the largest binary array are assigned to the xor result. (x^0 = x)
 */
long int *xor_binary_array(const long int *seq_bin1, const int seq_size1,
                           const long int *seq_bin2, const int seq_size2)
{

    // size of the binary array type used
    long int intsize = int_SIZE + 1;

    long int *s1, *s2;
    long int ss1, ss2;
    long int sbs1, sbs2;

    // Find the greater binary array, and rename them
    if (seq_size1 >= seq_size2)
    { // if s1 is greater or equal than s2
        s1 = seq_bin1;
        s2 = seq_bin2;
        ss1 = seq_size1 / intsize + ((seq_size1 / intsize) % intsize != 0);
        sbs1 = seq_size1;
        ss2 = seq_size2 / intsize + ((seq_size2 / intsize) % intsize != 0);
        sbs2 = seq_size2;
    }
    else
    { // else if s2 is greater or equal than s1
        s1 = seq_bin2;
        s2 = seq_bin1;
        ss1 = seq_size2 / intsize + ((seq_size2 / intsize) % intsize != 0);
        sbs1 = seq_size2;
        ss2 = seq_size1 / intsize + ((seq_size1 / intsize) % intsize != 0);
        sbs2 = seq_size1;
    }

    if (ss1 == 0)
        ss1 = 1;
    if (ss2 == 0)
        ss2 = 1;

    long int it = 0;

    // Allocate memory and verify it has been allocated
    long int * xor = NULL;
    xor = calloc(ss1, sizeof(*xor));
    if (!xor)
        return printf("ERROR: xor_binary_array: cannot allocate memory.\n"), NULL;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static, 64)
        for (it = 0; it < ss2 - 1; it++)
            xor[it] = s1[it] ^ s2[it];

#pragma omp master
        {
            // If one sequence is larger than the other, shift the last value of the smaller sequence toxor it with the other value.
            xor[ss2 - 1] = s1[ss2 - 1] ^ ((s2[ss2 - 1] << ((sbs1 - sbs2) % intsize)));
        }
        // Values from the largest binary array are assigned to the xor result. (x^0 = x)
#pragma omp for schedule(static, 64)
        for (it = ss2; it < ss1; it++)
            xor[it] = s1[it];
    }
    return xor;
}

/**
 * Popcount a binary array sequence.
 *
 * in : seq_bin : sequence in binary array format
 * in : seq_size : seq_bin length (number total of used bits)
 * out : bin_popcount : popcount of the seq : number of '1'
 *
 * Iterates on seq_bin and for each value, adds its popcount to bin_popcount.
 */
int popcount_binary_array(const long int *seq_bin, const long int seq_size)
{
    int bin_popcount = 0;

    // Find the size of the binary array
    long int array_size = seq_size / int_SIZE;
    if (seq_size % int_SIZE != 0)
        array_size++;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static, 64) reduction(+ \
                                               : bin_popcount)
        for (long i = 0; i < array_size; ++i)
            bin_popcount += __builtin_popcount(seq_bin[i]);
    }
    return bin_popcount;
}

/**
 * Retrieve a piece of the binary array sequence.
 *
 * in : seq_bin : sequence in binary array format
 * in : size : size of the binary requested
 * out : piece_seq_bin : the requested part of seq_bin, from pos_start and size.
 *
 * Iterates on seq_bin from pos_start, size times and gets for each iteration its binary value.
 */
long int *get_piece_binary_array(const long int *seq_bin, const unsigned long long pos_start, const unsigned long long pos_stop)
{
    // Find the size of the output
    long int array_size = (pos_stop - pos_start) / int_SIZE + ((pos_stop - pos_start) % int_SIZE != 0);

    // Allocate memory and verify it has been allocated
    long int *piece_seq_bin = NULL;
    piece_seq_bin = calloc(array_size, sizeof(*seq_bin));
    if (!piece_seq_bin)
        return printf("ERROR: get_piece_binary_array: cannot allocate memory.\n"), NULL;

        // stop position.

// Parse the binary array,
// from the bit at 'pos_start' position to 'pos_stop' position
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static, 64)
        for (int i = (int)pos_start; i < (int)pos_stop; i++)
        {
            int tmp = i - pos_start;
            change_binary_value(piece_seq_bin, (int)tmp, get_binary_value(seq_bin, (int)i));
        }
    }

    return piece_seq_bin;
}

/***************************************/
/******** DNA & GENES FUNCTION *********/
/***************************************/

//////////////// Convert to binary
/**
 * Convert a DNA base sequence to its binary array format.
 *
 * in : dna_seq : DNA sequence in char array format
 * in : size : dna_seq length = number of nucleotides (= number of letters)
 * out : seq : DNA sequence in binary array format
 *
 * Calls set_binary_array.
 */
long int *convert_to_binary(const char *dna_seq, size_t size)
{
    return set_binary_array(dna_seq, size);
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
char *binary_to_dna(long int *bin_dna_seq, const unsigned size)
{
    if (size % 2 != 0)
    {
        printf("Error: binary_to_aa : wrong binary size (%d). Must be odd.\nExit.\n", size);
        return NULL;
    }

    // Allocate memory and verify it has been allocated
    char *dna_seq = calloc((size / 2) + 1, sizeof(*dna_seq));
    if (!dna_seq)
        return printf("ERROR: binary_to_dna: cannot allocate memory.\n"), NULL;

    // Parse the binary array, two bits per iteration
    unsigned i;

    int nucl1, nucl2, index;
    char value;
#pragma omp parallel shared(bin_dna_seq, size, dna_seq, i) private(nucl1, nucl2, index, value)
    {
#pragma omp for schedule(static, 32)
        for (i = 0; i < size; i += 2)
        {
            // nucleotides = A, T, G, C
            nucl1 = get_binary_value(bin_dna_seq, i);
            nucl2 = get_binary_value(bin_dna_seq, i + 1);

            // get the ASCII value according to bits value
            value = bitstocharDNA[nucl2 + 2 * nucl1];

            index = i / 2;

            dna_seq[index] = value;
        }
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
char *generating_mRNA(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos)
{
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_mRNA: undefined sequence\n"), NULL;

    // Allocate memory and verify it has been allocated
    char *rna_seq = NULL;
    rna_seq = malloc(sizeof(*rna_seq) * ((stop_pos - start_pos) / 2) + 2);
    if (!rna_seq)
        return printf("ERROR: generating_mRNA: cannot allocate memory\n"), NULL;

    long int i;

    int nucl1, nucl2;
    char value;
    long index;

#pragma omp parallel shared(gene_seq, start_pos, stop_pos, i) private(nucl1, nucl2, value, index)
    {
        // Parse the binary DNA sequence, two bits per iteration
#pragma omp for schedule(static, 32)
        for (i = start_pos; i < stop_pos; i += 2)
        {

            // nucleotides = A, U, G, C
            nucl1 = get_binary_value(gene_seq, i);
            nucl2 = get_binary_value(gene_seq, i + 1);

            index = (i + 1 - start_pos) / 2;

            // get the ASCII value according to bits value
            value = bitstocharmRNA[nucl2 + 2 * nucl1];
            rna_seq[index] = value;
        }
    }
    // rna_seq[stop_pos / 2] = '\0';
    rna_seq[((stop_pos - start_pos) / 2) + 1] = '\0';
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
void detecting_genes(const long int *gene, const long int gene_size, gene_map_t *gene_map)
{
    gene_map->genes_counter = 0;

    // Check if memory ever have been allocated and allocate it if not
    if (!gene_map->gene_start || !gene_map->gene_end)
    {
        gene_map->gene_start = malloc(sizeof(*gene_map->gene_start) * MAX_GENES);
        gene_map->gene_end = malloc(sizeof(*gene_map->gene_end) * MAX_GENES);

        if (!gene_map->gene_start || !gene_map->gene_end)
        {
            printf("ERROR: detecting_genes: cannot allocate memory\n");
            return;
        }
    }

    // Start codon
    const __m128i AUG = _mm_set_epi16(0, 0, 1, 1, 0, 1, 0, 0);
    // Stop codons
    const __m128i UAA = _mm_set_epi16(1, 1, 0, 0, 0, 0, 0, 0);
    const __m128i UGA = _mm_set_epi16(1, 1, 0, 1, 0, 0, 0, 0);
    const __m128i UAG = _mm_set_epi16(1, 1, 0, 0, 0, 1, 0, 0);

    __m128i codon_to_test = _mm_setzero_si128();

    int start_pos = -1;

    long int i = 0;

    // Parse the binary array, and find all the start and stop codons
    while ((i + 6) <= gene_size)
    {
        // Each nucleotides can be A, U, G or C
        codon_to_test = _mm_set_epi16(get_binary_value(gene, i),
                                      get_binary_value(gene, i + 1),
                                      get_binary_value(gene, i + 2),
                                      get_binary_value(gene, i + 3),
                                      get_binary_value(gene, i + 4),
                                      get_binary_value(gene, i + 5),
                                      0, 0);

        const __m128i cmp_AUG = _mm_xor_si128(codon_to_test, AUG);

        const __m128i cmp_UAA = _mm_xor_si128(codon_to_test, UAA);
        const __m128i cmp_UAG = _mm_xor_si128(codon_to_test, UAG);
        const __m128i cmp_UGA = _mm_xor_si128(codon_to_test, UGA);

        // If a start pos and a stop pos doesn't exist, search for AUG
        if (_mm_testz_si128(cmp_AUG, cmp_AUG))
        {
            // if AUG, it's the start of a gene
            start_pos = i;
            i += 6;
        }
        else if ((start_pos != -1) && (_mm_testz_si128(cmp_UAA, cmp_UAA) || _mm_testz_si128(cmp_UAG, cmp_UAG) || _mm_testz_si128(cmp_UGA, cmp_UGA)))
        {
            // It's the end of a gene
            // If a start pos and an stop pos has been found, a gene exists so we save it in the struc
            gene_map->gene_start[gene_map->genes_counter] = start_pos;
            gene_map->gene_end[gene_map->genes_counter] = i + 5;

            gene_map->genes_counter++;

            start_pos = -1;
            i += 6;
        }
        else
            i += 2;
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
char *generating_amino_acid_chain(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos)
{
    long int codon_size = 6;
    // Check the input argument
    if (!gene_seq)
        return printf("ERROR: generating_amino_acid_chain: undefined sequence\n"), NULL;
    if ((stop_pos - start_pos) % 3 != 0)
        return NULL;

    // Allocate memory and verify it has been allocated
    char *aa_seq = calloc(sizeof(char), sizeof(*aa_seq) * ((stop_pos - start_pos) / codon_size) + 1);
    if (!aa_seq)
        return printf("ERROR: generating_amino_acid_chain: cannot allocate memory\n"), NULL;

    unsigned temp = 0;

    long int k;
    int tmp, pow_bit, get_bin;

    long int i, size = stop_pos - start_pos;

#pragma omp parallel shared(gene_seq, start_pos, codon_size, size, aa_seq, i) private(k, temp, tmp, pow_bit, get_bin)
    {
        // Parse the binary array, six bits by six (to parse three nucleotides per three)
#pragma omp for schedule(static, 64 / codon_size)
        for (i = start_pos; i < stop_pos; i += codon_size)
        {

            temp = (i - start_pos) / codon_size;

            // Get the decimal value of the 6 bits
            tmp = 0;
            pow_bit = 5;
            for (k = i; k < i + codon_size; k++)
            {
                get_bin = get_binary_value(gene_seq, k);
                tmp += get_bin << pow_bit;
                pow_bit--;
            }

            // Get the corresponding protein from the lookup table
            aa_seq[temp] = LUT[tmp];
        }
    }
    aa_seq[size / codon_size] = '\0';
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
void detecting_mutations(const long int *gene_seq, const unsigned long long start_pos, const unsigned long long stop_pos,
                         mutation_map mut_m)
{
    unsigned long long detect_mut = 0;
    unsigned long long size_sequence = stop_pos - start_pos; // Counting size of GC sequence
    unsigned short tmp_start_mut = 0;                        // stock start mutation
    unsigned cmp = 0;                                        // counter of all mutation zones
    const __m128i G = _mm_set_epi64x((long)0, (long)1);
    const __m128i C = _mm_set_epi64x((long)1, (long)0);

    // long size = start_pos + size_sequence;
    // Parse the binary array, from the 'start_pos' bit to the end
    for (unsigned long long i = start_pos; i < stop_pos; i += 2)
    {

        // each nucleotides can be  A, U, G or C
        const long bit1 = (long)get_binary_value(gene_seq, i);
        const long bit2 = (long)get_binary_value(gene_seq, i + 1);

        const __m128i nucl = _mm_set_epi64x(bit1, bit2);
        const __m128i cmp_G = _mm_xor_si128(nucl, G);
        const __m128i cmp_C = _mm_xor_si128(nucl, C);

        // Increment detect_mut if find a C or G nucl
        if (_mm_testz_si128(cmp_G, cmp_G) || _mm_testz_si128(cmp_C, cmp_C))
        {
            if (detect_mut == 0)
                tmp_start_mut = i - start_pos;
            detect_mut += 2;
        }
        // Put detect_mut to 0 if find a A or T nucl
        else
        {
            // Check if previous GC sequence is a probable mutation zone
            if (detect_mut >= (size_sequence / 5))
            {
                mut_m.start_mut[cmp] = tmp_start_mut;

                mut_m.end_mut[cmp] = (i)-start_pos;

                mut_m.size[cmp] = detect_mut - 1;

                cmp++;
            }
            detect_mut = 0;
        }
    }
    // Check if ending sequence is a probable mutation zone
    if (detect_mut >= (size_sequence / 5))
    {
        mut_m.start_mut[cmp] = tmp_start_mut;
        mut_m.end_mut[cmp] = size_sequence;
        mut_m.size[cmp] = detect_mut - 1;

        cmp++;
    }
}

#ifdef __AVX512VPOPCNTDQ__
__attribute__((target("avx512vpopcntdq", "avx512f"))) float calculating_matching_score(long int *seq1, const int seq_size1,
                                                                                       long int *seq2, const int seq_size2)
{
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    long int *big_seq = NULL, *small_seq = NULL;
    int len_big = 0, len_small = 0, diff_len = 0;

    // Find and rename the biggest/smallest sequences
    // !! size = number of bits, length = size of the array (= size/int_SIZE)
    if (seq_size1 >= seq_size2)
    {
        big_seq = seq1;
        small_seq = seq2;

        len_big = seq_size1 / int_SIZE;
        if (seq_size1 % int_SIZE != 0)
            len_big++;

        len_small = seq_size2 / int_SIZE;
        if (seq_size2 % int_SIZE != 0)
            len_small++;
    }
    else
    {
        big_seq = seq2;
        small_seq = seq1;

        len_big = seq_size2 / int_SIZE;
        if (seq_size2 % int_SIZE != 0)
            len_big++;

        len_small = seq_size1 / int_SIZE;
        if (seq_size1 % int_SIZE != 0)
            len_small++;
    }

    // Parse both seqences from the end
    /**
     * (Need to add some "0" at the beginning, if too small)
     *
     * big =   l1    l2    l3     l4     l5
     * small =  0     0    l3'    l4'    l5'
     * ---------------------------------------------
     * xor =   l1    l2   l3^l3' l4^l4' l5^l5'
     */

    int pop = 0;
    int j = len_small;
    for (int i = len_big; i > 0; i -= 8)
    {
        const long x0 = j >= 0 ? (long)small_seq[j - 1] : (long)0;
        const long x1 = j - 1 >= 0 ? (long)small_seq[j - 2] : (long)0;
        const long x2 = j - 2 >= 0 ? (long)small_seq[j - 3] : (long)0;
        const long x3 = j - 3 >= 0 ? (long)small_seq[j - 4] : (long)0;
        const long x4 = j - 4 >= 0 ? (long)small_seq[j - 5] : (long)0;
        const long x5 = j - 5 >= 0 ? (long)small_seq[j - 6] : (long)0;
        const long x6 = j - 6 >= 0 ? (long)small_seq[j - 7] : (long)0;
        const long x7 = j - 7 >= 0 ? (long)small_seq[j - 8] : (long)0;

        const long y0 = j >= 0 ? (long)big_seq[i - 1] : (long)0;
        const long y1 = j - 1 >= 0 ? (long)big_seq[i - 2] : (long)0;
        const long y2 = j - 2 >= 0 ? (long)big_seq[i - 3] : (long)0;
        const long y3 = j - 3 >= 0 ? (long)big_seq[i - 4] : (long)0;
        const long y4 = j - 4 >= 0 ? (long)big_seq[i - 5] : (long)0;
        const long y5 = j - 5 >= 0 ? (long)big_seq[i - 6] : (long)0;
        const long y6 = j - 6 >= 0 ? (long)big_seq[i - 7] : (long)0;
        const long y7 = j - 7 >= 0 ? (long)big_seq[i - 8] : (long)0;

        const __m512i xor = _mm512_xor_epi64(_mm512_set_epi64(y0, y1, y2, y3, y4, y5, y6, y7),
                                             _mm512_set_epi64(x0, x1, x2, x3, x4, x5, x6, x7));

        const __m512i popcnt = _mm512_popcnt_epi64(xor);

        pop += _mm512_reduce_add_epi64(popcnt);
        j -= 8;
    }

    const float denom = 1.0 / (float)(len_big * int_SIZE);
    const float F = pop ? (float)pop * 100.0 * denom : 0.0;

    return 100.0 - F;
}
#else
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
float calculating_matching_score(long int *seq1, const int seq_size1,
                                 long int *seq2, const int seq_size2)
{
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    // First step: apply the xor operation between both arrays
    long int * xor = NULL;
    xor = xor_binary_array(seq1, seq_size1, seq2, seq_size2);
    // xor_size = max size between 'seq_size1' and 'seq_size2'
    int xor_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;

    // Second step: count the number of bits whose value is 1 on the result
    int pop = popcount_binary_array(xor, xor_size);

    // Last step: compute the percentage
    float y = ((float)pop * 100.0) / (float)xor_size;
    return 100.0 - y;
}
#endif

#define MAX(x, y) ((x) > (y) ? (x) : (y))

//////////////// Print Similarity Matrix
/**
 * Prints the similarity matrix
 *
 * in : F : The similarity matrix in COL MAJOR. Its dimensions are `m*n`.
 * in : A : Horizontal char array. Its size is `m-1`
 * in : B : Vertical char array. Its size is `n-1`
 * in : m : Number of columns of `F`. Also the `A` length + 1
 * in : n : Number of lines of `F`. Also the `B` length + 1
 * in : k : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 * in : l : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 */
void print_sim_mat(int* F, char A [], char B [], int m, int n, int k, int l) {
    printf("%3s  %3s  ", "*", "-");
    for (int j = 0; j < m; j++)
        printf("%3c  ", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%3c  ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 3d  \033[0m", F[i * m + j]);
            else
                printf("% 3d  ", F[i * m + j]);
        }
        printf("\n");
    }
}

//////////////// Calculate Similarity Matrix (anti-diagonal)
/**
 * Calculates the similarity matrix of two sequences.
 * Iteration method is anti-diagonal
 * Sequences are stored as char arrays.
 * One extra space for the similarity matrix is required, its returned dimensions are (`A` length+1)*(`B` length+1`
 *
 * in : A : Horizontal char array.
 * in : B : Vertical char array.
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * out : F : The similarity matrix in COL MAJOR. <!> Its dimensions are (`A` length+1)*(`B` length+1).
 */
int* calculate_scoring_matrix_antidiag(char A [], char B [], int match, int mismatch, int gap) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < m; i++)
        F[i] = gap * i;
    for (int j = 0; j < n; j++)
        F[j * m] = gap * j;

    int max_l, i, j;
    for (int k = 1; k < m + n - 2; k++) {
        if (k > (m - 1)) {
            max_l = n - (k - (m - 1)) - 1;
            if (max_l > m - 1) max_l = m - 1;
        }
        else max_l = k;

#pragma omp parallel for schedule(guided)
        for (int l = 1; l <= max_l; l++) {
            if (k > (m - 1)) {
                i = k - (m - 1) + l;
                j = (m - 1) - l + 1;
            }
            else {
                i = l;
                j = k - i + 1;
            }
            diag = F[(i - 1) * m + j - 1] + (A[j - 1] == B[i - 1] ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    }
    return F;
}

//////////////// String Insert
/**
 * Insert a char at beginning of a char array.
 *
 * in/out : dest   : Destination string, to which we insert the `source` char
 * in     : source : Char to be inerted at beginning of `dest`.
 */
void string_insert(char** dest, char* source) {
    if (*dest == NULL) {
        *dest = calloc(sizeof(char), 2);
        memcpy(*dest, source, 1);
        return;
    }
    char* temp = calloc(sizeof(char), strlen(*dest) + 2);
    memcpy(temp, source, 1);
    strcat(temp, *dest);
    free(*dest);
    *dest = temp;
}

//////////////// Align
/**
 * Backtrace the similarity matrix, and prints the two sequences aligned.
 * Returns the max similarity score of the two sequences.
 *
 * in  : F : The similarity matrix in COL MAJOR. Its dimensions are given by (`A` length + 1)*(`B` length + 1).
 * in  : A : Horizontal char array.
 * in  : B : Vertical char array.
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * in  : print : Enables (1) or Disables (0) alignment output. Returns max similarity score either ways.
 * out : Max similarity score of the `A` and `B` sequences.
 */
int align(int* F, char A [], char B [], int match, int mismatch, int gap, int print) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    if (print) {
        int i = m - 1;
        int j = n - 1;

        char* A_aligned = NULL;
        char* B_aligned = NULL;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && F[i * m + j] == F[(i - 1) * m + j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch)) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, B + (j - 1));
                i--;
                j--;
            }
            else if (i > 0 && F[i * m + j] == F[(i - 1) * m + j] + gap) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, "-");
                i--;
            }
            else {
                // else if (j > 0 && F[i * m + j] == F[i * m + j - 1] + gap) {
                string_insert(&A_aligned, "-");
                string_insert(&B_aligned, B + (j - 1));
                j--;
            }
        }
        printf("%s\n", A_aligned);
        printf("%s\n", B_aligned);
    }
    return F[m * n - 1];
}

//////////////// Needleman-Wunsch (anti-diag)
/**
 * Basic Needleman-Wunsch application.
 * Given two char arrays, returns the maximum similarity score.
 *
 * in  : A : Horizontal char array.
 * in  : B : Vertical char array.
 * in  : m : Length of `A`.
 * in  : n : Length of `B`.
 * out : Max similarity score of the `A` and `B` sequences.
 */
int needleman_wunsch_antidiag(char A[], char B[]) {
    int m = strlen(A);
    int n = strlen(B);

    int match = 1;
    int mismatch = -1;
    int gap = -1;

    int* F = calculate_scoring_matrix_antidiag(A, B, match, mismatch, gap);

    // int score = align(F, A, B, match, mismatch, gap, 1);
    // printf("Max score : %d\n", score);
    int score = F[(m + 1) * (n + 1) - 1];
    free(F);
    return score;
}

/* Count files in a directory
 *
 */
int countfiles()
{
    int count = 0;
    struct dirent *entry;

    DIR *dir = opendir("./fastas/");

    while ((entry = readdir(dir)) != NULL)
    {
        if (strstr(entry->d_name, ".fasta"))
            count++;
    }
    closedir(dir);
    return count;
}

void insert_list(node_t **head, long int *data, int size)
{
    if (data == NULL)
        printf("Error data \n");

    node_t *tmp_node = (node_t *)malloc(sizeof(node_t));
    node_t *last = *head;

    tmp_node->seq = malloc(sizeof(long) * (size));

    for (int i = 0; i < size; i++)
    {
        tmp_node->seq[i] = data[i];
    }

    tmp_node->size = size;
    tmp_node->next = NULL;

    if (*head == NULL)
    {
        tmp_node->prev = NULL;
        *head = tmp_node;
    }
    else
    {
        while (last->next != NULL)
            last = last->next;

        last->next = tmp_node;
        tmp_node->prev = last;
    }
}

int readfiles(int comm_size, int nbseq){
    int align_rank = comm_size - 1;

    int nb = countfiles();
    char **content = malloc(sizeof(char *) * nb);
    FILE *fp;

    if (output == 1){
        struct stat st = {0};
        if (stat("./output/", &st) == -1)
            mkdir("./output/", 0700);

        fp = fopen("./output/rapport_bin.html", "w+");
        if (fp == NULL){
            printf("Cannot open file \n");
            exit(0);
        }
        fprintf(fp, "<html>\n<head><style>\n th, td {\n        font - size : 10px; \n}\n.title {\n        font - size : 15px; \n}\ntable, th, td {\n    border:\n        1px solid black;\n        border - collapse : collapse;\n        border - style : dashed;\n}\n.title {\n        border - style : dashed dashed dashed solid;\n        padding - left : 1 %% ;\n}\ntable {\n    width:\n        90 %% ;\n        margin - left : 5 %% ;\n}\n\n\ndetails > summary {\n    padding:\n        4px;\n    width:\n        200px;\n        background - color : #eeeeee;\n    border:\n        none;\n        box - shadow : 1px 1px 2px #bbbbbb;\n    cursor:\n        help;\n}\n</style>\n</head>\n");
    }
    DIR *dir;
    FILE *input;
    struct dirent *file;

    char folder_name[50] = "./fastas/";

    // Open the directory which contain all the fastas files
    if ((dir = opendir(folder_name)) == NULL)
        return printf("Error: Can't open fastas folder\n"), -1;

    int i = 0;
    node_t *head = NULL;

    // Iterate if a file exists in this directory
    while (((file = readdir(dir)) != NULL) && (i < nbseq)){
        // Skip directories (linux)
        if (file->d_type == DT_DIR)
            continue;
        if ((!strcmp(file->d_name, ".")) && (!strcmp(file->d_name, "..")))
            continue;

        printf("Reading file : %s\n", file->d_name);

        // Get filepath
        char *filepath = NULL;
        filepath = malloc(sizeof(char) * (strlen(folder_name) + strlen(file->d_name)));
        strcpy(filepath, folder_name);
        strcat(filepath, file->d_name);

        // Get file size (will be used for content allocation)
        struct stat st;
        stat(filepath, &st);
        long filesize = st.st_size;

        content[i] = (char *)malloc(filesize * sizeof(char));

        // Open fasta file
        if ((input = fopen(filepath, "r")) == NULL)
            return printf("Error: Can't open fastas file %s\n", filepath), -1;

        char *line;
        size_t len = 0;
        ssize_t read;

        // Skipping first line, as it's only fasta's metadata.
        read = getline(&line, &len, input);

        // Init content[i] with second line (first data line)
        content[i] = malloc((filesize - read - 1) * sizeof(char));
        getline(&line, &len, input);
        line[strcspn(line, "\n") - 1] = '\0';
        strcpy(content[i], line);

        // Concat each lines in content[i], while toggling newline.
        while ((read = getline(&line, &len, input)) != -1){
            line[strcspn(line, "\n") - 1] = '\0';
            strcat(content[i], line);
        }
        free(line);
        fclose(input);

        int recv = i % (comm_size - 2) + 1;

        if (output)
            fprintf(fp, "<details><summary> %s </summary>\n<a href=\"sequences/rank_%d_%d_bin.html\"> %s </a></details>\n", file->d_name, recv, i / comm_size, file->d_name);

        // printf("%d) sending to %d\n", rank, recv);
        MPI_Send(content[i], strlen(content[i]), MPI_CHAR, recv, 0, MPI_COMM_WORLD);

        // Sends seq for alignment purposes.
        MPI_Send(content[i], strlen(content[i]), MPI_CHAR, align_rank, 4, MPI_COMM_WORLD);
        i++;
    }

    if (closedir(dir) == -1)
        return printf("Error close dir\n"), -1;

    MPI_Send(&i, 1, MPI_INT, align_rank, 5, MPI_COMM_WORLD);

    i = 0;

    // Sends the sequence count to receiver. Ends communication from MASTER to receiver
    for (int j = 1; j < comm_size - 1; j++)
        MPI_Send(&i, 1, MPI_INT, j, 1, MPI_COMM_WORLD);

    int cont = 1;

    MPI_Status status;

    // printf("%d) comm_size = %d", rank, comm_size);

    // while (cont < comm_size) {
    do
    {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == 3)
        {
            MPI_Recv(&i, 1, MPI_INT, status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cont++;
            // printf("cont = %d comm_size = %d \n", cont, comm_size);
        }
        else
        {
            // printf("%d receive from %d | tag : %d\n", rank, status.MPI_SOURCE, status.MPI_TAG);

            int count;

            MPI_Get_count(&status, MPI_LONG, &count);
            long int *tmp = malloc(sizeof(long int) * count);

            MPI_Recv(tmp, count, MPI_LONG, status.MPI_SOURCE, 2, MPI_COMM_WORLD, &status);

            insert_list(&head, tmp, count);
            free(tmp);
        }
    } while (cont < comm_size - 1);

    printf("Rank 0 ended receiving data\n");

    node_t *seq1;
    node_t *seq2;

    seq1 = head;
    // printf("%ld\n", seq1->seq[0]);

    FILE *comp;

    if (output)
    {

        fprintf(fp, "<details><summary> %s </summary>\n<a href=\"comp/comp_bin.html\"> Matching score </a></details>\n", "");

        struct stat st = {0};
        if (stat("./output/comp", &st) == -1)
        {
            mkdir("./output/comp", 0700);
        }

        comp = fopen("./output/comp_bin.html", "w+");

        if (fp == NULL)
        {
            printf("Cannot open file \n");
            exit(0);
        }
        fprintf(comp, "<html>\n<head><style>\n th, td {\n        font - size : 10px; \n}\n.title {\n        font - size : 15px; \n}\ntable, th, td {\n    border:\n        1px solid black;\n        border - collapse : collapse;\n        border - style : dashed;\n}\n.title {\n        border - style : dashed dashed dashed solid;\n        padding - left : 1 %% ;\n}\ntable {\n    width:\n        90 %% ;\n        margin - left : 5 %% ;\n}\n\n\ndetails > summary {\n    padding:\n        4px;\n    width:\n        200px;\n        background - color : #eeeeee;\n    border:\n        none;\n        box - shadow : 1px 1px 2px #bbbbbb;\n    cursor:\n        help;\n}\n</style>\n</head>\n");

        fprintf(comp, "<table>\n<tbody>\n<tr>\n<td class = \"title\">Sequence 1</td>\n<td class = \"title\">Sequence 2</td>\n<td class = \"title\">Score</td></tr>\n");
    }

    while (seq1 != NULL)
    {
        seq2 = seq1->next;
        while (seq2 != NULL)
        {
            // printf("%ld\n", seq2->seq[0]);
            float calc = calculating_matching_score(seq1->seq, seq1->size * 32, seq2->seq, seq2->size * 32);
            if (output)
                fprintf(comp, "<tr><td>%s</td><td>%s</td><td>%f</td></tr>", binary_to_dna(seq1->seq, seq1->size * 32), binary_to_dna(seq2->seq, seq2->size * 32), calc);
            // printf("Matching score : %f\n", calc);
            seq2 = seq2->next;
        }

        seq1 = seq1->next;
    }

    // Free everything
    free(content);

    if (output)
    {
        fprintf(comp, "</table> </body></html>");
        fclose(comp);

        fprintf(fp, "</html>");
        fclose(fp);
    }
    printf("Rank 0 terminated\n");
    return 0;
}

void alignment_work(int rank) {
    int real_nb_seqs = 0;
    MPI_Status sta;

    char** seq = malloc(sizeof(char*) * 200);

    int cont = 1;
    int i = 0;
    int flag;
    while (cont) {
        MPI_Status status;

        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        if (flag) {
            if (status.MPI_TAG == 5) {
                MPI_Status sta;
                MPI_Recv(&real_nb_seqs, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &sta);
                if (i == real_nb_seqs)
                    cont = 0;
            }
            else if (status.MPI_TAG == 4) {
                int count = 0;
                MPI_Get_count(&status, MPI_CHAR, &count);
                seq[i] = (char*)malloc(sizeof(char) * count);
                MPI_Recv(seq[i], count, MPI_CHAR, 0, 4, MPI_COMM_WORLD, &sta);
                i++;

                if (real_nb_seqs && i == real_nb_seqs)
                    cont = 0;
            }
        }
    }

    for (int i = 0; i < real_nb_seqs - 1; i++) {
        needleman_wunsch_antidiag(seq[i], seq[i+1]);
        // printf("*%d\tAlignment score (%d:%d) : %d\n", rank, i, i+1, score);
        printf("Rank %d - %.0f%% (%d/%d)\n", rank, (float)(i + 1) / (real_nb_seqs-1) * 100, (i + 1), real_nb_seqs-1);
        free(seq[i]);
    }
    free(seq[real_nb_seqs-1]);
    free(seq);
}

void process_work(int rank)
{
    int flag;
    int cont = 1;

    int nb = countfiles();
    char *seq[nb];

    int i = 0;
    int count = 0;

    while (cont)
    {
        MPI_Status status;

        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        if (flag)
        {
            if (status.MPI_TAG == 1)
            {
                MPI_Status sta;
                MPI_Recv(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &sta);
                cont = 0;
            }
            else
            {
                MPI_Get_count(&status, MPI_CHAR, &count);
                seq[i] = (char *)malloc(sizeof(char) * count);
                MPI_Recv(seq[i], count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
                i++;
            }
        }
    }

    for (int j = 0; j < i; j++)
    {
        gene_map_t gene_map;
        long *seq_bin;
        long len_seq;

        seq_bin = convert_to_binary(seq[j], strlen(seq[j]));

        len_seq = strlen(seq[j]) * 2;

        gene_map.gene_start = malloc(sizeof(*gene_map.gene_start) * strlen(seq[j]) * int_SIZE);
        gene_map.gene_end = malloc(sizeof(*gene_map.gene_end) * strlen(seq[j]) * int_SIZE);

        detecting_genes(seq_bin, len_seq, &gene_map);

        // printf("*%d\tGenes detected : %lld\n", rank, gene_map.genes_counter);
        FILE *fp;

        if (output){
            struct stat st = {0};

            if (stat("./output/sequences", &st) == -1)
                mkdir("./output/sequences", 0700);

            char *rank_c, *nb_c;
            char name_f[50];

            asprintf(&rank_c, "%d", rank);
            asprintf(&nb_c, "%d", j);

            strcat(strcpy(name_f, "./output/sequences/rank_"), rank_c);
            strcat(strcat(name_f, "_"), nb_c);
            strcat(name_f, "_bin.html");

            fp = fopen(name_f, "w+");

            if (fp == NULL){
                printf("Cannot open file \n");
                exit(0);
            }
            fprintf(fp, "<html>\n<head><style>\n th, td {\n        font - size : 10px; \n}\n.title {\n        font - size : 15px; \n}\ntable, th, td {\n    border:\n        1px solid black;\n        border - collapse : collapse;\n        border - style : dashed;\n}\n.title {\n        border - style : dashed dashed dashed solid;\n        padding - left : 1 %% ;\n}\ntable {\n    width:\n        90 %% ;\n        margin - left : 5 %% ;\n}\n\n\ndetails > summary {\n    padding:\n        4px;\n    width:\n        200px;\n        background - color : #eeeeee;\n    border:\n        none;\n        box - shadow : 1px 1px 2px #bbbbbb;\n    cursor:\n        help;\n}\n</style>\n</head>\n");
        }
        for (unsigned long long k = 0; k < gene_map.genes_counter; k++)
        {
            mutation_map mut_m;

            long *genes = get_piece_binary_array(seq_bin, gene_map.gene_start[k], gene_map.gene_end[k]);

            char *amino = generating_amino_acid_chain(seq_bin, gene_map.gene_start[k], gene_map.gene_end[k]);
            // if (amino != NULL)
            //     printf("*%d\tamino acid chain = %s\n", rank, amino);

            mut_m.size = malloc(sizeof(*mut_m.size) * ((gene_map.gene_end[k] - gene_map.gene_start[k]) / 5) * int_SIZE);
            mut_m.start_mut = malloc(sizeof(*mut_m.start_mut) * ((gene_map.gene_end[k] - gene_map.gene_start[k]) / 5) * int_SIZE);
            mut_m.end_mut = malloc(sizeof(*mut_m.end_mut) * ((gene_map.gene_end[k] - gene_map.gene_start[k]) / 5) * int_SIZE);

            char *mrna = generating_mRNA(seq_bin, gene_map.gene_start[k], gene_map.gene_end[k]);
            // printf("*%d\tMRNA = %s\n", rank, mrna);

            detecting_mutations(seq_bin, gene_map.gene_start[k], gene_map.gene_end[k], mut_m);

            int size_m = (gene_map.gene_end[k] - gene_map.gene_start[k]) / int_SIZE;
            if (!size_m)
                size_m = 1;

            if (output)
            {
                fprintf(fp, "<table>\n<tbody>\n<tr>\n<td class = \"title\">Sequence</td>\n<td class = \"title\">MRNA</td>\n<td class = \"title\">Chain</td>\n<td class = \"title\">Mutation</td>\n</tr>\n");

                fprintf(fp, "<tr><td>%s</td>", binary_to_dna(genes, (gene_map.gene_end[k] - gene_map.gene_start[k]) + 3));
                fprintf(fp, "<td>%s</td>\n", mrna);
                if (amino != NULL)
                    fprintf(fp, "<td>%s</td>\n", amino);
                else
                    fprintf(fp, "<td>NONE</td>\n");

                if (mut_m.size[0] == 0)
                    fprintf(fp, "<td>NONE</td>\n");
                else
                {
                    fprintf(fp, "<td>");
                    int c = 0;
                    while ((mut_m.size[0] != 0) && (c < 5))
                    {
                        fprintf(fp, "[%lu,%lu]", mut_m.start_mut[c], mut_m.end_mut[c]);
                        c++;
                    }
                    fprintf(fp, "</td>\n");
                }

                fprintf(fp, "</tbody>\n</table>\n");
            }

            // printf("*%d\tSending %d genes to 0\n", rank, gene_map.genes_counter);
            MPI_Send(genes, size_m, MPI_LONG, 0, 2, MPI_COMM_WORLD);

            free(mut_m.end_mut);
            free(mut_m.size);
            free(mut_m.start_mut);
        }

        if (output){
            fprintf(fp, "</html>");
            fclose(fp);
        }
        printf("Rank %d - %.0f%% (%d/%d)\n", rank, (float)(j + 1) / i*100, (j + 1), i);
    }

    MPI_Send(&cont, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
}

void launch(int affichage, int sequences)
{
    output = affichage;
    int RANK_MASTER = 0;

    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(NULL, NULL);

    //   int rank;
    int comm_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (rank == RANK_MASTER)
    {
        readfiles(comm_size, sequences);
    }
    else if (rank == comm_size - 1)
        alignment_work(rank);
    else
    {
        process_work(rank);
    }

    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
}