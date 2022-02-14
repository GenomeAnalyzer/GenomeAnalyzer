#include <stdint.h>
#include <string.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "../headers/gene_bin.h"
#include "../src/gene_bin.c"

// static void test_get_binary_value(void ** state){
//     // Test if the algorithm is OK
//     // 1 = 0000000000000000000000000000001
//     assert_int_equal(1, get_binary_value((long int []){1}, 0));
//     //9350764 = 001101100111010101110001
//     assert_int_equal(0, get_binary_value((long int []){1}, 17));
//     // Pour chaque binaire de 00000 à 11111, vérifier chaque bit.
//     for(int i =0; i<32;i++)
//     {
//         assert_int_equal(i%2, get_binary_value((long int  []){i}, 0));
//         assert_int_equal(i/2%2, get_binary_value((long int  []){i}, 1));
//         assert_int_equal(i/4%2, get_binary_value((long int  []){i}, 2));
//         assert_int_equal(i/8%2, get_binary_value((long int  []){i}, 3));
//         assert_int_equal(i/16%2, get_binary_value((long int  []){i}, 4));
//     }
// }

static void test_convert_to_binary(void** state) {
    // Test aa to binary conversions
    // --- Test all valid letters
    // A 00 T 11 C 10 G 01

    // Check all valid letters.
    // Binary sequence is inverted.
    long int* array = calloc(1, sizeof(long int));
    convert_to_binary(array, "ATCG\n", 5);
    assert_int_equal(0x3900000000000000, array[0]);

    array = calloc(1, sizeof(long int));
    convert_to_binary(array, "ATCGATCGATCG\n", 13);
    assert_int_equal(0x3939390000000000, array[0]);

    // Give a different size than the char sequence.
    array = calloc(1, sizeof(long int));
    convert_to_binary(array, "ATCG\n", 5);
    assert_int_equal(0x3900000000000000, array[0]);

    array = calloc(1, sizeof(long int));
    convert_to_binary(array, "ATCG\n", 6);
    assert_int_equal(0x39, array[0]); // Expect 00 for low order bit

    char* seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCG\nATCGATCGATCGATCGATCGATCGATCGATCG\n";
    array = calloc(3, sizeof(long int));
    convert_to_binary(array, seq_char, 33);
    assert_int_equal(0x3939393939393939, array[0]);

    array = calloc(3, sizeof(long int));
    convert_to_binary(array, seq_char, 66);
    assert_int_equal(0x3939393939393939, array[1]);

    // Check of a random sequence
    seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCT\nTCGAGACCTTCGA\n";
    unsigned seq_size = 47;
    array = calloc(2, sizeof(long int));
    convert_to_binary(array, seq_char, seq_size);

    long int seq_sol[2] = {5402369836413063467, (long unsigned int)59845604 << 38};
    assert_int_equal(seq_sol[0], array[0]);
    assert_int_equal(seq_sol[1], array[1]);

    // Test whether the function correctly detects errors:
    // --- Unknown letter in sequence
    convert_to_binary(array, "AJ\n", 10);
    assert_int_equal(-1, array[0]);

    free(array);
}

// static void test_binary_to_dna(void ** state){
//     // Test binary to aa conversions

//     // --- Test all conversion
//     assert_string_equal("ATCG", binary_to_dna((long int []){156}, 8));

//     char* seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCG";
//     long int seq_size = 32;
//     long int* seq_bin;
//     char* seq_test = NULL;
//     seq_test = calloc(seq_size, sizeof(char));
//     char* seq_new = NULL;
//     int* ptr;

//     for(long int i = 1; i < seq_size; i++)
//     {
//         seq_bin = convert_to_binary(seq_char, i);
//         seq_new = binary_to_dna(seq_bin, 2*i);
//         ptr = (int*) realloc(seq_test, sizeof(char)*i);
//         // Check the realloc worked.
//         assert_ptr_not_equal(NULL, ptr);
//         // Copy the sequence in the test sequence, to keep only the sequence needed for the test
//         memcpy(seq_test, seq_char, i);
//         assert_string_equal(seq_test, seq_new);
//     }

//     // Test whether the function correctly detects errors:
//     // --- Wrong size bin_dna_seq
//     assert_ptr_equal(NULL, binary_to_dna((long int []){0}, 3));
// }

// static void test_generating_mRNA(void ** state){
//     // Test if the algorithm is OK
//     //9350764 = 001101100111010101110001
//     char* seq_char = NULL;
//     assert_string_equal("AUGCGUGGGUAG", generating_mRNA((long int []){9350764}, 0, 24));

//     seq_char = generating_mRNA((long int []) { 1821290092, 18263 }, 0, 128);
//     assert_string_equal("AUGCGUGGGUAGAUGC", generating_mRNA((long int []) { 1821290092 }, 0, 32));
//     assert_string_equal("AUGCGUGGGUAGAUGCAAAAAAAAAAAAAAAA", generating_mRNA((long int []) { 1821290092 }, 0, 64));
//     assert_string_equal("AUGCGUGGGUAGAUGCAAAAAAAAAAAAAAAAGUGGGUAGAAAAAAAAAAAAAAAAAAAAAAAA", generating_mRNA((long int []) { 1821290092, 18263  }, 0, 128));

//     seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
//     char* expected_char = "AUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG";

//     long int seq_size = 64;
//     long int* seq_bin;
//     char* seq_mRNA = NULL;

//     seq_bin = convert_to_binary(seq_char, 2 * seq_size);
//     seq_mRNA = generating_mRNA(seq_bin, 0, 2*seq_size);
//     assert_string_equal(expected_char, seq_mRNA);

//     expected_char = "AUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG";
//     seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCG";
//     seq_size = 32;
//     char* seq_test = NULL;
//     seq_test = calloc(seq_size, sizeof(char));
//     char* seq_new = NULL;
//     int* ptr;

//     for (long int i = 1; i < seq_size; i++)
//     {
//         seq_bin = convert_to_binary(seq_char, i);
//         seq_new = generating_mRNA(seq_bin, 0, 2 * i);
//         ptr = (int*)realloc(seq_test, sizeof(char) * i);
//         // Check the realloc worked.
//         assert_ptr_not_equal(NULL, ptr);
//         // Copy the sequence in the test sequence, to keep only the sequence needed for the test
//         memcpy(seq_test, expected_char, i);
//         assert_string_equal(seq_test, seq_new);
//     }

//     // Test whether the function correctly detects errors:
//     // --- NULL error
//     assert_ptr_equal(NULL, generating_mRNA(NULL, 0, 0));
// }


static void test_detecting_genes(void ** state){
    gene_map_t gene_map;
    gene_map.gene_start = malloc(sizeof(*gene_map.gene_start) * 3);
    gene_map.gene_end = malloc(sizeof(*gene_map.gene_end) * 3);
    gene_map.gene_start[0] = 0;
    gene_map.gene_end[0] = 0;
    gene_map.gene_start[1] = 0;
    gene_map.gene_end[1] = 0;
    gene_map.gene_start[2] = 0;
    gene_map.gene_end[2] = 0;

    // Test if the algorithm is OK in a basic case: xxxAUGxxxxUAAxxx
    //                                              AGC AUG AGGGCC UAA CGU
    // The algorithm should detect one gene from the start codon to the stop codon

    long int* seq_bin = calloc(1, sizeof(long int));
    convert_to_binary(seq_bin, "AGCATGAGGGCCTAACGT\n", 19);
    // printf("seq_bin : %s \n", seq_bin);
    detecting_genes(seq_bin, 18*2, &gene_map);
    assert_int_equal(6, gene_map.gene_start[0]);
    assert_int_equal(29, gene_map.gene_end[0]);
    assert_int_equal(1, gene_map.genes_counter);

    // Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
    //                                                         AGC AUG AGGGCC AUG CGAACG UAA CGU
    // The algorithm should detect one gene from the 2nd start codon to the stop codon
    seq_bin = calloc(1, sizeof(long int));
    convert_to_binary(seq_bin, "AGCATGAGGGCCATGCGAACGTAACGT\n", 28);
    detecting_genes(seq_bin, 27*2, &gene_map);
    assert_int_equal(24, gene_map.gene_start[0]);
    assert_int_equal(47, gene_map.gene_end[0]);
    assert_int_equal(1, gene_map.genes_counter);

    // Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
    //                                                         AGC AUG CACGCG UAA GCACTG UAA CGU
    // The algorithm should detect one gene from the start codon to the first stop codon
    seq_bin = calloc(1, sizeof(long int));
    convert_to_binary(seq_bin, "AGCATGCACGCGTAAGCACTGTAACGT\n", 28);
    detecting_genes(seq_bin, 27*2, &gene_map);
    assert_int_equal(6, gene_map.gene_start[0]);
    assert_int_equal(29, gene_map.gene_end[0]);
    assert_int_equal(1, gene_map.genes_counter);

    // Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
    //                                                                     CGCCGCGCCGCGGGCG
    // The algorithm should not detect any genes
    gene_map.gene_start[0] = 0;
    gene_map.gene_end[0] = 0;
    seq_bin = calloc(1, sizeof(long int));
    convert_to_binary(seq_bin, "CGCCGCGCCGCGGGCG\n", 17);
    detecting_genes(seq_bin, 16*2, &gene_map);
    assert_int_equal(0, gene_map.gene_start[0]);
    assert_int_equal(0, gene_map.gene_end[0]);
    assert_int_equal(0, gene_map.genes_counter);

    // Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
    //                                                      AGC AUG GCGCAC UAG CGCCCG AUG CUGGGG UAA CGU
    // The algorithm should detect two genes
    seq_bin = calloc(1, sizeof(long int));
    convert_to_binary(seq_bin, "AGCATGGCGCACTAGCGCCCGATGCTGGGGTA\nACGT\n", 38);
    detecting_genes(seq_bin, 36*2, &gene_map);
    assert_int_equal(6, gene_map.gene_start[0]);
    assert_int_equal(29, gene_map.gene_end[0]);
    assert_int_equal(42, gene_map.gene_start[1]);
    assert_int_equal(65, gene_map.gene_end[1]);
    assert_int_equal(2, gene_map.genes_counter);

    free(gene_map.gene_start);
    free(gene_map.gene_end);
    free(seq_bin);
}

// static void test_generating_aa_chain(void ** state){
//     // Test if the algorithm is OK

//     // AA -> SEQ -> BINARY -> INVERT BINARY -> CALL FUNCTION
//     // K      K      N      N      R      R      S      S      T      T     
//     // AAA    AAG    AAC    AAT    AGA    AGG    AGC    AGT    ACA    ACG         
//     // 000000 000001 000010 000011 000100 000101 000110 000111 001000 001001   
//     assert_string_equal("KKNNRRSSTT", generating_amino_acid_chain((long int []) {0b100100000100111000011000101000001000110000010000100000000000}, 0, 60));
//     // T      T      I      M      I      I      E      E      D      D
//     // ACC    ACT    ATA    ATG    ATC    ATT    GAA    GAG    GAC    GAT
//     // 001010 001011 001100 001101 001110 001111 010000 010001 010010 010011
//     assert_string_equal("TTIMIIEEDD", generating_amino_acid_chain((long int []) { 0b110010010010100010000010111100011100101100001100110100010100 }, 0, 60));
//     // G      G      A      A      A      A      V      V      V      V 
//     // GGC    GGT    GCA    GCG    GCC    GCT    GTA    GTG    GTC    GTT         
//     // 010110 010111 011000 011001 011010 011011 011100 011101 011110 011111   
//     assert_string_equal("GGAAAAVVVV", generating_amino_acid_chain((long int []) { 0b111110011110101110001110110110010110100110000110111010011010 }, 0, 60));
//     // Q      Q      H      H      R      R      R      R      P      P
//     // CAA    CAG    CAC    CAT    CGA    CGG    CGC    CGT    CCA    CCG
//     // 100000 100001 100010 100011 100100 100101 100110 100111 101000 101001
//     assert_string_equal("QQHHRRRRPP", generating_amino_acid_chain((long int []) { 0b100101000101111001011001101001001001110001010001100001000001 }, 0, 60));
//     // L      L      L      L      O      O      Y      Y      O      W     
//     // CTA    CTG    CTC    CTT    TAA    TAG    TAC    TAT    TGA    TGG    
//     // 101100 101101 101110 101111 110000 110001 110010 110011 110100 110101
//     assert_string_equal("LLLLOOYYOW", generating_amino_acid_chain((long int []) { 0b101011001011110011010011100011000011111101011101101101001101 }, 0, 60));
//     // C      C      S      S      S      S      L      L      F      F
//     // TGC    TGT    TCA    TCG    TCC    TCT    TTA    TTG    TTC    TTT
//     // 110110 110111 111000 111001 111010 111011 111100 111101 111110 111111
//     assert_string_equal("CCSSSSLLFF", generating_amino_acid_chain((long int []) { 0b111111011111101111001111110111010111100111000111111011011011 }, 0, 60));
//     // G      G      P      P
//     // GGA    GGG    CCC    CCT
//     // 010100 010101 101010 101011
//     assert_string_equal("GGPP", generating_amino_acid_chain((long int []) { 0b110101010101101010001010 }, 0, 24));


//     //  --- Test all the amino acid
//     // (alphabetic order of the above sequences.)
//     long int* seq_bin = convert_to_binary("AAAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTGAAGAGGACGATGGCGGTGCAGCGGCCGCTGTAGTGGTCGTTCAACAGCACCATCGACGGCGCCGTCCACCGCTACTGCTCCTTTAATAGTACTATTGATGGTGCTGTTCATCGTCCTCTTTATTGTTCTTTGGAGGGCCCCCT", 384);
//     char* aa_chain = NULL;
//     aa_chain = generating_amino_acid_chain(seq_bin, 0, 384);
//     assert_string_equal("KKNNRRSSTTTTIMIIEEDDGGAAAAVVVVQQHHRRRRPPLLLLOOYYOWCCSSSSLLFFGGPP", aa_chain);

//     // Test whether the function correctly detects errors:
//     // --- NULL error
//     assert_ptr_equal(NULL, generating_amino_acid_chain(NULL, 0, 0));
//     // --- invalid value in gene_seq
//     // assert_ptr_not_equal("MRGO",
//     //                      generating_amino_acid_chain((long int []){1821290092, 18263}, 24));
// }


// static void test_detecting_mutations(void ** state){
//     mutation_map M;
//     unsigned short nb_mutations = 6;
//     M.size = calloc(nb_mutations, sizeof(unsigned long));
//     M.start_mut = calloc(nb_mutations, sizeof(unsigned long));
//     M.end_mut = calloc(nb_mutations, sizeof(unsigned long));

//     //A : 00
//     //T : 11
//     //C : 10
//     //G : 01

//     //GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
//     long int* seq_bin = convert_to_binary("GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG", 30);
//     //Test if sequence 10 to 23 is a mutation zone and no other mutation zone
//     M.size[1]=0;
//     M.start_mut[1]=0;
//     M.end_mut[1]=0;
//     detecting_mutations(seq_bin, 0, 60, M);
//     // detecting_mutations((long int []){261725162, 97523700}, 0, 60, M);
//     //First mutation is updated with right values
//     assert_int_equal(13,M.size[0]);
//     assert_int_equal(10,M.start_mut[0]);
//     assert_int_equal(24,M.end_mut[0]);
//     //No other mutations, should not be updated
//     assert_int_equal(0,M.size[1]);
//     assert_int_equal(0,M.start_mut[1]);
//     assert_int_equal(0,M.end_mut[1]);

//     //GTTTTGCAAACGTTAAAGGTTTGAAAGGTG = {261102590, 97523700}
//     seq_bin = convert_to_binary("GTTTTGCAAACGTTAAAGGTTTGAAAGGTG", 30);
//     //Test if no mutation in this sequence
//     M.size[0]=0;
//     M.start_mut[0]=0;
//     M.end_mut[0]=0;
//     detecting_mutations((long int []){261102590, 97523700}, 0, 60, M);
//     //No possible mutation zones detected, should not be updated
//     assert_int_equal(0,M.size[0]);
//     assert_int_equal(0,M.start_mut[0]);
//     assert_int_equal(0,M.end_mut[0]);

//     //GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
//     seq_bin = convert_to_binary("GGGCCGTTCCGCCCATAGGCCCGGCTAAGA", 30);
//     //Test with 3 mutation zones in this sequence
//     M.size[3]=0;
//     M.start_mut[3]=0;
//     M.end_mut[3]=0;
//     detecting_mutations(seq_bin, 0, 60, M);
//     //First mutation is updated with right values
//     assert_int_equal(11,M.size[0]);
//     assert_int_equal(0,M.start_mut[0]);
//     assert_int_equal(12,M.end_mut[0]);
//     //Second mutation is updated with right values
//     assert_int_equal(11,M.size[1]);
//     assert_int_equal(16,M.start_mut[1]);
//     assert_int_equal(28,M.end_mut[1]);
//     //Third mutation is updated with right values
//     assert_int_equal(15,M.size[2]);
//     assert_int_equal(34,M.start_mut[2]);
//     assert_int_equal(50,M.end_mut[2]);
//     //No other mutations, should not be updated
//     assert_int_equal(0,M.size[3]);
//     assert_int_equal(0,M.start_mut[3]);
//     assert_int_equal(0,M.end_mut[3]);


//     // Test three mutations in sequence
//     seq_bin = convert_to_binary("GGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTAAAGGTG", 147);
//     detecting_mutations(seq_bin, 0, 147*2, M);

//     for(int i = 0; i < 3; i++){
//         assert_int_equal(85, M.size[i]);
//         assert_int_equal(10 + 90*i, M.start_mut[i]);
//         assert_int_equal(95 + 90 * i+1, M.end_mut[i]);
//     }

//     free(M.size);
//     free(M.start_mut);
//     free(M.end_mut);
// }

static void test_calculating_matching_score(void ** state){
    // Test if the algorithm is OK
    // --- With same size
    //  GACCCGAC = 0100101010010010 = {18770} -> 3 popcounts
    //  GGCCAGGC = 0101101000010110 = {26714}
    assert_float_equal(100.0-3*100.0/16.0,
                    calculating_matching_score(
                      (long int []){18770}, 0, 16,
                      (long int []){26714}, 0, 16),
                    0);
    // --- With different size
    //  GACCCGAC =   0100101010010010     = {18770}
    //  TTTCAGGCTC = 11111110000101101110 = {485503}
    assert_float_equal(100.0-9*100.0/20.0,
                    calculating_matching_score(
                      (long int []){18770}, 0, 16,
                      (long int []){485503}, 0, 20),
                    0);
    //  TTTCAGGCTT = 11111110000101101111 = {1009791} -> 15 popcounts
    //  GACCTTCGA =  010010101111100100   = {40786}
    assert_float_equal(100.0-15*100.0/20,
                    calculating_matching_score(
                      (long int []){1009791}, 0, 20,
                      (long int []){40786}, 0, 18),
                    0);
    assert_float_equal(100.0-15*100.0/20,
                    calculating_matching_score(
                      (long int []){1009791}, 0, 20,
                      (long int []){40786}, 0, 16),
                    0);

    // --- With more than 1 element
    //  GACCCGACGACCCGACGACCCGACGACCCGACTTTCAGGCTT -> 27 popcounts
    //  GGCCAGGCGGCCAGGCGGCCAGGCGGCCAGGCGACCTTCGA
    assert_float_equal(100.0-27*100.0/84.0,
                    calculating_matching_score(
                      (long int []){5283365930625288000, 1009791}, 0, 84,
                      (long int []){7519437265355565000, 40786}, 0, 82),
                    0);

    // Test whether the function correctly detects errors:
    // --- NULL error
    assert_float_equal(-1.0, calculating_matching_score(NULL, 0, 0, NULL, 0, 0), 0);
}

int main(void) {
    int result = 0;
    const struct CMUnitTest tests[] = {
        // BINARIES ARRAYS FUNCTIONS
        // cmocka_unit_test(test_get_binary_value),
        // DNA & GENES FUNCTIONS
        cmocka_unit_test(test_convert_to_binary),
        // cmocka_unit_test(test_binary_to_dna),
        // cmocka_unit_test(test_generating_mRNA),
        cmocka_unit_test(test_detecting_genes),
        // cmocka_unit_test(test_generating_aa_chain),
        // cmocka_unit_test(test_detecting_mutations),
        cmocka_unit_test(test_calculating_matching_score),
    };
    result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

    return result;
}
