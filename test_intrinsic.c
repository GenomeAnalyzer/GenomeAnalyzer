#include <nmmintrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "versions/bin/headers/gene_bin.h"
#include "rdtsc/rdtsc.h"

#define MAX 36
#define MAX_LINE 150
#define MAX_ADN 100000

/*
void load_gene(char *filename, char *seq_char)
{
    char temp[MAX_LINE] = "";
    FILE *file;
    file = fopen(filename, "r");
    if(!file)
        printf("ERROR: load_gene: cannot open file %s\n", filename);

    if(file)
    {
        fgets(temp, MAX_LINE, file); 
        while(fgets(temp, MAX_LINE, file) != NULL)
        {
            temp[strlen(temp)-1] = '\0';
            strcat(seq_char, temp);
        }
    }
    fclose(file);
}
*/

unsigned long long load_from_file(char *filename, char *seq_char)
{
	FILE *file;
	file = fopen(filename, "r");
	if(!file)
		printf("ERROR: load_from_file: cannot open file 'seq.txt'\n");

	unsigned long long i = 0;
	while((seq_char[i] = fgetc(file)) != EOF)
		i++;
	printf("size seq_char  = %llu\n", i-1);

	fclose(file);
	return i-1;
}

void detecting_genes_2(const long int *gene, const long int gene_size, gene_map_t* gene_map) {
        gene_map->genes_counter = 0;

    // Check if memory ever have been allocated and allocate it if not
    if(!gene_map->gene_start || !gene_map->gene_end){
        gene_map->gene_start = malloc(sizeof(*gene_map->gene_start) * MAX_GENES);
        gene_map->gene_end = malloc(sizeof(*gene_map->gene_end) * MAX_GENES);

        if (!gene_map->gene_start || !gene_map->gene_end){
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
   //Parse the binary array, and find all the start and stop codons
    while ((i + 6) <= gene_size) {
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

        //If a start pos and a stop pos doesn't exist, search for AUG
        if (_mm_testz_si128(cmp_AUG, cmp_AUG)) {
        //if AUG, it's the start of a gene
            start_pos = i;
            i += 6;
        }
        else if ((start_pos != -1 ) 
	         && (_mm_testz_si128(cmp_UAA, cmp_UAA)
	             || _mm_testz_si128(cmp_UAG, cmp_UAG)
                     || _mm_testz_si128(cmp_UGA, cmp_UGA))) {
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

void test_m128i()
{
    int test = 0;
    __m128i v0 = _mm_setzero_si128();
    __m128i v1 = _mm_setzero_si128();
    __m128i vcmp = _mm_setzero_si128();

    v0 = _mm_set_epi16(0,1,0,1,0,1,1,1);
    v1 = _mm_set_epi16(1,1,1,1,1,1,0,0);

    vcmp = _mm_xor_si128(v0, v1);
    test = _mm_testz_si128(vcmp, vcmp);

    if(test)
        printf("equal: %d\n\n", test);
    else
        printf("different: %d\n\n", test);
}

unsigned long long test1_genes(unsigned long long seq_char_size, long int *seq_long)
{
    gene_map_t g;
    g.gene_start = malloc(sizeof(*g.gene_start) * seq_char_size * int_SIZE);
    g.gene_end   = malloc(sizeof(*g.gene_end) * seq_char_size * int_SIZE);

    detecting_genes(seq_long, 2 * seq_char_size, &g);
    //printf("%llu genes\n\n", g.genes_counter);

    free(g.gene_start);
    free(g.gene_end);
    return g.genes_counter;
}

unsigned long long test2_genes(unsigned long long seq_char_size, long int *seq_long)
{
    gene_map_t g;
    g.gene_start = malloc(sizeof(*g.gene_start) * seq_char_size * int_SIZE);
    g.gene_end   = malloc(sizeof(*g.gene_end) * seq_char_size * int_SIZE);

    detecting_genes_2(seq_long, 2 * seq_char_size, &g);
    //printf("%llu genes\n\n", g.genes_counter);

    free(g.gene_start);
    free(g.gene_end);
    return g.genes_counter;
}

void detecting_mutations_2(const long int *gene_seq, const long int start_pos,
                           const long int size_sequence, mutation_map mut_m) {
    long int detect_mut = 0;  //Counting size of GC sequence
    unsigned short tmp_start_mut = 0;   //stock start mutation
    unsigned cmp = 0;   //counter of all mutation zones

    const __m128i G = _mm_set_epi64x((long) 0, (long) 1);
    const __m128i C = _mm_set_epi64x((long) 1, (long) 0);

    long size = start_pos + size_sequence;
    //Parse the binary array, from the 'start_pos' bit to the end
    for (long int i = start_pos; i < size; i += 2) {

        // each nucleotides can be  A, U, G or C
        const long bit1 = (long) get_binary_value(gene_seq, i);
        const long bit2 = (long) get_binary_value(gene_seq, i + 1);

        const __m128i nucl = _mm_set_epi64x(bit1, bit2);
        const __m128i cmp_G = _mm_xor_si128(nucl, G);
        const __m128i cmp_C = _mm_xor_si128(nucl, G);

        //Increment detect_mut if find a C or G nucl
        if (_mm_testz_si128(cmp_G, cmp_G) || _mm_testz_si128(cmp_C, cmp_C)) {
            if(detect_mut == 0)
                tmp_start_mut = i-start_pos;
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

        //TMP//TMP//TMP
        cmp++;
        printf("cmp = %d\n", cmp);
        //TMP//TMP//TMP
    }
    else
        printf("cmp = %d\n", cmp);
}

unsigned long long test1_mutations(unsigned long long seq_char_size, long int *seq_long)
{
    mutation_map mut_m;
    mut_m.size =      malloc(sizeof(*mut_m.size) * (seq_char_size/5) * int_SIZE);
    mut_m.start_mut = malloc(sizeof(*mut_m.start_mut) * (seq_char_size/5) * int_SIZE);
    mut_m.end_mut =   malloc(sizeof(*mut_m.end_mut) * (seq_char_size/5) * int_SIZE);

    detecting_mutations(seq_long, 0, 2 * seq_char_size, mut_m);
    //printf("%llu genes\n\n", g.genes_counter);

    free(mut_m.start_mut);
    free(mut_m.end_mut);
    return 0;
}

unsigned long long test2_mutations(unsigned long long seq_char_size, long int *seq_long)
{
    mutation_map mut_m;
    mut_m.size =      malloc(sizeof(*mut_m.size) * (seq_char_size/5) * int_SIZE);
    mut_m.start_mut = malloc(sizeof(*mut_m.start_mut) * (seq_char_size/5 + 1) * int_SIZE);
    mut_m.end_mut =   malloc(sizeof(*mut_m.end_mut) * (seq_char_size/5 + 1) * int_SIZE);

    detecting_mutations_2(seq_long, 0, 2 * seq_char_size, mut_m);
    //printf("%llu genes\n\n", g.genes_counter);

    free(mut_m.start_mut);
    free(mut_m.end_mut);
    return 0;
}

/* ORIGINEL (seulement ajouté un printf au milieu) */
float calculating_matching_score_0(const long int *seq1, long int start_pos1,const int seq_size1,
                                 const long int *seq2, long int start_pos2,const int seq_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    // First step: apply the xor operation between both arrays
    long int *seq1tmp;
    seq1tmp = get_piece_binary_array(seq1, start_pos1, seq_size1);

    long int *seq2tmp;
    seq2tmp = get_piece_binary_array(seq2, start_pos2, seq_size2);

    long int *xor = NULL;
    xor = xor_binary_array(seq1tmp, seq_size1, seq2tmp, seq_size2);
    // xor_size = max size between 'seq_size1' and 'seq_size2'
    int xor_size = seq_size1 >= seq_size2 ? seq_size1 : seq_size2;

    //Second step: count the number of bits whose value is 1 on the result
    int pop = popcount_binary_array(xor, xor_size);

    printf("pop = %d\n", pop);

    //Last step: compute the percentage
    float y = ((float)pop * 100.0) / (float)xor_size;
    return 100.0 - y;
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
__attribute__((target("avx512vpopcntdq", "avx512f")))
float calculating_matching_score_2(long int *seq1,const int seq_size1,
                                   long int *seq2,const int seq_size2) {
    // Check the input argument
    if (!seq1 || !seq2)
        return printf("ERROR: calculating_matching_score: undefined sequence\n"), -1.0;

    long int *big_seq = NULL, *small_seq = NULL;
    int len_big = 0, len_small = 0, diff_len = 0;

    // Find and rename the biggest/smallest sequences
    // !! size = number of bits, length = size of the array (= size/int_SIZE)
    if(seq_size1 >= seq_size2){
        big_seq = seq1;
        small_seq = seq2;

        len_big = seq_size1 / int_SIZE;
        if (seq_size1 % int_SIZE != 0) len_big++;

        len_small = seq_size2 / int_SIZE;
        if (seq_size2 % int_SIZE != 0) len_small++;
    } 
    else{
        big_seq = seq2;
        small_seq = seq1;

        len_big = seq_size2 / int_SIZE;
        if (seq_size2 % int_SIZE != 0) len_big++;

        len_small = seq_size1 / int_SIZE;
        if (seq_size1 % int_SIZE != 0) len_small++;
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
    for (int i = len_big; i > 0; i-= 8)
    {
        const long x0 = j   >=0 ? (long)small_seq[j - 1] : (long)0;
        const long x1 = j-1 >=0 ? (long)small_seq[j - 2] : (long)0;
        const long x2 = j-2 >=0 ? (long)small_seq[j - 3] : (long)0;
        const long x3 = j-3 >=0 ? (long)small_seq[j - 4] : (long)0;
        const long x4 = j-4 >=0 ? (long)small_seq[j - 5] : (long)0;
        const long x5 = j-5 >=0 ? (long)small_seq[j - 6] : (long)0;
        const long x6 = j-6 >=0 ? (long)small_seq[j - 7] : (long)0;
        const long x7 = j-7 >=0 ? (long)small_seq[j - 8] : (long)0;

        const long y0 = j   >=0 ? (long)big_seq[i - 1] : (long)0;
        const long y1 = j-1 >=0 ? (long)big_seq[i - 2] : (long)0;
        const long y2 = j-2 >=0 ? (long)big_seq[i - 3] : (long)0;
        const long y3 = j-3 >=0 ? (long)big_seq[i - 4] : (long)0;
        const long y4 = j-4 >=0 ? (long)big_seq[i - 5] : (long)0;
        const long y5 = j-5 >=0 ? (long)big_seq[i - 6] : (long)0;
        const long y6 = j-6 >=0 ? (long)big_seq[i - 7] : (long)0;
        const long y7 = j-7 >=0 ? (long)big_seq[i - 8] : (long)0;

        const __m512i xor = _mm512_xor_epi64(_mm512_set_epi64(y0, y1, y2, y3, y4, y5, y6, y7),
                                             _mm512_set_epi64(x0, x1, x2, x3, x4, x5, x6, x7));

        const __m512i popcnt = _mm512_popcnt_epi64(xor);

        pop += _mm512_reduce_add_epi64(popcnt);
        j -= 8;
    }

    printf("pop = %d\n", pop);
 
    const float denom = 1.0 / (float)(len_big * int_SIZE);
    const float F = pop ? (float)pop * 100.0 * denom : 0.0;

    return 100.0 - F;

}

int main(int argc, char *argv[])
{
    // ------------------- //
    // ---- TEST AVX2 ---- //
    // ------------------- //
    test_m128i();


    // ------------------------------- //
    // ---- TESTS DETECTING GENES ---- //
    // ------------------------------- //
    char *seq_char = NULL, *seq_char2 = NULL;
    unsigned long long seq_char_size = 0, seq_char_size2 = 0;
    long int *seq_long = NULL, *seq_long2 = NULL;

    seq_char = malloc(sizeof(*seq_char) * MAX_ADN);
    if(!seq_char)
	    return printf("error malloc\n"), 1;
    seq_char_size = load_from_file("seq_test_intrinsic_2.txt", seq_char);
    seq_long = convert_to_binary(seq_char, seq_char_size);

    seq_char2 = malloc(sizeof(*seq_char2) * MAX_ADN);
    if(!seq_char2)
        return printf("error malloc\n"), 1;
    seq_char_size2 = load_from_file("seq_test_intrinsic.txt", seq_char2);
    seq_long2 = convert_to_binary(seq_char2, seq_char_size2);
 
    double cycles = 0.0;
    double boucles = 0.0;
    unsigned long long cnt = 0;

    // real detecting genes
    while(cycles <= 500)
    {
    	double before = (double)rdtsc();
        cnt = test1_genes(seq_char_size, seq_long);
    	double after = (double)rdtsc();

    	cycles += (after - before);
    	boucles++;
    }
    printf("real detecting genes : %lf cycles, %0.1lf boucles, genes cnt = %llu\n", cycles / boucles, boucles, cnt);

    // intrinsic detecting genes
    cycles = 0.0;
    boucles = 0.0;
    cnt = 0;
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        cnt = test2_genes(seq_char_size, seq_long);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("intrinsic detecting genes : %lf cycles, %0.1lf boucles, genes cnt = %llu\n\n", cycles / boucles, boucles, cnt);



    // ----------------------------------- //
    // ---- TESTS DETECTING MUTATIONS ---- //
    // ----------------------------------- //

    cycles = 0.0;
    boucles = 0.0;
    cnt = 0;

    // real detecting mutations
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        cnt = test1_mutations(seq_char_size2, seq_long2);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("real detecting mutations : %lf cycles, %0.1lf boucles, mutations cnt = %llu\n", cycles / boucles, boucles, cnt);

    // intrinsic detecting genes
    cycles = 0.0;
    boucles = 0.0;
    cnt = 0;
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        cnt = test2_mutations(seq_char_size2, seq_long2);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("intrinsic detecting mutations : %lf cycles, %0.1lf boucles, mutations cnt = %llu\n\n", cycles / boucles, boucles, cnt);

    // ------------------------------------ //
    // ------- TESTS MATCHING SCORE ------- //
    // ------------------------------------ //

    cycles = 0.0;
    boucles = 0.0;
    float MS = 0.0;

    // real calculating_matching_score
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        MS = calculating_matching_score_0(seq_long, 0, 2*seq_char_size, seq_long2, 0, 2*seq_char_size2);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("real calculating_m_s : %lf cycles, %0.1lf boucles, matching score = %f\n", cycles / boucles, boucles, MS);


    cycles = 0.0;
    boucles = 0.0;
    MS = 0.0;
    // new calculating_matching_score
    while(cycles <= 500)
    {
        double before = (double)rdtsc();
        MS = calculating_matching_score_2(seq_long, 2*seq_char_size, seq_long2, 2*seq_char_size2);
        double after = (double)rdtsc();

        cycles += (after - before);
        boucles++;
    }
    printf("new calculating_m_s : %lf cycles, %0.1lf boucles, matching score = %f\n", cycles / boucles, boucles, MS);

    free(seq_char);
    free(seq_long);

    free(seq_char2);
    free(seq_long2);

    return 0;
}
