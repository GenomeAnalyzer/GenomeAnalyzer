#include <stdint.h>
#include <string.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "../headers/gene_bin.h"
#include "../src/gene_bin.c"

static void test_get_binary_value(void** state) {

  // 1 = 0000000000000000000000000000001
  assert_int_equal(1, get_binary_value((long int []) { 1 }, 0));

  // Requesting a bit out of the binary range should returns 0
  assert_int_equal(0, get_binary_value((long int []) { 0b11111 }, 6));

  // For every bit from 0b000000000 to 0b111111111, verify each bit, for each value.
  // We loop on all the 9bits numbers (2^0 -> 2^8-1)
  // The value of the bit at position x (from the right) corresponds to i / 2^x % 2
  for (int value = 0; value < 128;i++) {
    assert_int_equal(value % 2, get_binary_value((long int []) { value }, 0));
    assert_int_equal(value / 2 % 2, get_binary_value((long int []) { value }, 1));
    assert_int_equal(value / 4 % 2, get_binary_value((long int []) { value }, 2));
    assert_int_equal(value / 8 % 2, get_binary_value((long int []) { value }, 3));
    assert_int_equal(value / 16 % 2, get_binary_value((long int []) { value }, 4));
    assert_int_equal(value / 32 % 2, get_binary_value((long int []) { value }, 5));
    assert_int_equal(value / 64 % 2, get_binary_value((long int []) { value }, 6));
  }

  // Requesting bits outside of the range should returns 0 (31 bits for a long int.)
  assert_int_equal(1, get_binary_value((long int []) { 0b1111111111111111111111111111111111111111 }, 31));
  assert_int_equal(0, get_binary_value((long int []) { 0b1111111111111111111111111111111111111111 }, 32));
}

static void test_change_binary_value(void ** state){
  // Test if the algorithm is OK

  long int seq_bin = 0b1111111111111111111111111111111111111111111111111111111111111111;

  for(int i = 0; i < 64; i++) change_binary_value(&seq_bin, i, 0);
  // printf("seq_bin : %ld\n", seq_bin);
  assert_int_equal(seq_bin, 0);
  for(int i = 0; i < 64; i++) change_binary_value(&seq_bin, i, 1);
  assert_int_equal(seq_bin, 0b1111111111111111111111111111111111111111111111111111111111111111);

  // 85 = 1010101
  seq_bin = 85;
  // Invert values of seq_bin
  for (int i = 0;i < 7;i++) change_binary_value(&seq_bin, i, i%2);
  // 42 = 0101010
  assert_int_equal(seq_bin, 42);
  // Invert again values of seq_bin.
  for (int i = 0;i < 7;i++) change_binary_value(&seq_bin, i, (i+1)%2);
  assert_int_equal(seq_bin, 85);
}

static void test_convert_to_binary(void** state) {
  // Test aa to binary conversions
  // --- Test all valid letters
  // A 00 T 11 C 10 G 01

  // Check all valid letters.
  // Binary sequence is inverted.
  assert_int_equal(0b10011100, convert_to_binary("ATCG", 4)[0]);
  assert_int_equal(0b100111001001110010011100, convert_to_binary("ATCGATCGATCG", 12)[0]);

  // Give a different size than the char sequence.
  assert_int_equal(0b011100, convert_to_binary("ATCG", 3)[0]);
  assert_int_equal(0b0010011100, convert_to_binary("ATCG", 5)[0]); // Expect 00 for high order bit

  char* seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
  assert_int_equal(0b10011100100111001001110010011100, convert_to_binary(seq_char, 16)[0]);
  assert_int_equal(0b1001110010011100100111001001110010011100100111001001110010011100, convert_to_binary(seq_char, 32)[0]);

  // Check for multiple positions in array
  assert_int_equal(0b1001110010011100100111001001110010011100100111001001110010011100, convert_to_binary(seq_char, 64)[0]);

  // Check of a random sequence
  seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
  unsigned seq_size = 45;

  long int* seq_bin = NULL;
  seq_bin = convert_to_binary(seq_char, seq_size);

  long int seq_sol[2] = { -3131702537379864750, -9223372036849555181 };
  assert_int_equal(seq_sol[0], seq_bin[0]);
  assert_int_equal(seq_sol[1], seq_bin[1]);

  // Test whether the function correctly detects errors:
  // --- Unknown letter in sequence
  long int* res2 = convert_to_binary("AK", 10);
  assert_ptr_equal(NULL, res2[2]);
}

static void test_set_binary_array(void ** state){
  // Test if the algorithm is OK

  assert_int_equal(0b111100111111, set_binary_array("TTTATT", 6)[0]);
  assert_int_equal(0b11111111, set_binary_array("TTTT", 4)[0]);
  assert_int_equal(0b0, set_binary_array("AAAA", 4)[0]);
  assert_int_equal(0b01010101, set_binary_array("CCCC", 4)[0]);
  assert_int_equal(0b01101100, set_binary_array("ATGC", 4)[0]);
  assert_int_equal(0b00111001, set_binary_array("CGTA", 4)[0]);

  char* seq_char = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

  long int  pow = 1;
  long int  pow2 = 1;
  for (unsigned i = 1; i < 64; i++) {
    if (i <= 32) {
      pow *= 4;
      // printf("bin[0] => %ld, %ld\n", pow - 1, set_binary_array(seq_char, i)[0]);
      assert_int_equal(pow-1, set_binary_array(seq_char, i)[0]);
    }
    else {
      // printf("bin[1] => %ld, %ld\n", 0b1000000000000000000000000000000000000000000000000000000000000000 + pow2*2 - 1, set_binary_array(seq_char, i)[1]);
      assert_int_equal(0b1000000000000000000000000000000000000000000000000000000000000000 + pow2 * 2 - 1, set_binary_array(seq_char, i)[1]);
      pow2 *= 4;
    }
  }
}

static void test_xor_binary_array(void ** state){
  // Test if the algorithm is OK
  long int *xor = NULL;
  xor = xor_binary_array((long int  []) {42}, 7, (long int  []) {85}, 7);
  assert_int_equal(xor[0], 127);

  char *seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
  unsigned seq_size = 45;
  long int  *seq_bin = NULL;
  seq_bin = set_binary_array(seq_char, seq_size);

  char *seq_char2 = "GACCTTTTTTTTTTTTTCTTCGA";
  unsigned seq_size2 = 23;
  long int  *seq_bin2 = NULL;
  seq_bin2 = set_binary_array(seq_char2, seq_size2);

  seq_char = "TTTAATTTAATTTAATTTAATTTAATTTAATTT";
  seq_char = "TATAAATATTATAAATAT";
  seq_size = 18;
  seq_bin = set_binary_array(seq_char, seq_size);
  // printf("=> seq_bin = %d\n", seq_bin[0]);
  // 0b11 11 11 11 11 11
  xor = xor_binary_array(seq_bin, 2 * seq_size, seq_bin, 2 * seq_size);
  // printf("=> xor = %d\n", xor [0]);

  // // long int  *xor = NULL;
  // xor = xor_binary_array(seq_bin, 2 * seq_size, seq_bin2, 2 * seq_size2);
  // // xor_size = max size
  // int xor_size = seq_size >= seq_size2 ? seq_size : seq_size2;

  // long int  xor_sol[3] = {2101911378, 364800657, 2934};

  // for (int i = 0; i < 3; ++i)
  //   assert_int_equal(xor[i], xor_sol[i]);

  free(seq_bin);
  free(seq_bin2);
  free(xor);
}

static void test_popcount_binary_array(void ** state){
  // Test if the algorithm is OK
      // 1 = 0000000000000000000000000000001
  assert_int_equal(1, popcount_binary_array((long int []){1}, 31));
      //9350764 = 001101100111010101110001
  assert_int_equal(13, popcount_binary_array((long int []){9350764}, 24));

  // Test for all binaries from 00000 to 11111
  int popc_expected_result = 0, popc_result = 0;
  for (int i = 0; i < 32; i++){
    popc_expected_result = i%2 + i/2%2 + i/4%2 + i/8%2 + i/16%2;
    popc_result = popcount_binary_array((long int  []) { i }, 31);
    assert_int_equal(popc_expected_result, popc_result);
  }
}

static void test_binary_to_dna(void ** state){
  // Test binary to aa conversions

  // --- Test all conversion
  assert_string_equal("ATCG", binary_to_dna((long int []){156}, 8));

  char* seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCG";
  long int seq_size = 32;
  long int* seq_bin;
  char* seq_test = NULL;
  seq_test = calloc(seq_size, sizeof(char));
  char* seq_new = NULL;
  int* ptr;

  for(long int i = 1; i < seq_size; i++){
    seq_bin = convert_to_binary(seq_char, i);
    seq_new = binary_to_dna(seq_bin, 2*i);
    ptr = (int*) realloc(seq_test, sizeof(char)*i);
    // Check the realloc worked.
    assert_ptr_not_equal(NULL, ptr);
    // Copy the sequence in the test sequence, to keep only the sequence needed for the test
    memcpy(seq_test, seq_char, i);
    assert_string_equal(seq_test, seq_new);
  }

  // Test whether the function correctly detects errors:
  // --- Wrong size bin_dna_seq
  assert_ptr_equal(NULL, binary_to_dna((long int []){0}, 3));
}

static void test_generating_mRNA(void ** state){
  // Test if the algorithm is OK
      //9350764 = 001101100111010101110001
  char* seq_char = NULL;
  assert_string_equal("AUGCGUGGGUAG", generating_mRNA((long int []){9350764}, 0, 24));

  seq_char = generating_mRNA((long int []) { 1821290092, 18263 }, 0, 128);
  assert_string_equal("AUGCGUGGGUAGAUGC", generating_mRNA((long int []) { 1821290092 }, 0, 32));
  assert_string_equal("AUGCGUGGGUAGAUGCAAAAAAAAAAAAAAAA", generating_mRNA((long int []) { 1821290092 }, 0, 64));
  assert_string_equal("AUGCGUGGGUAGAUGCAAAAAAAAAAAAAAAAGUGGGUAGAAAAAAAAAAAAAAAAAAAAAAAA", generating_mRNA((long int []) { 1821290092, 18263  }, 0, 128));

  seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
  char* expected_char = "AUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG";

  long int seq_size = 64;
  long int* seq_bin;
  char* seq_mRNA = NULL;

  seq_bin = convert_to_binary(seq_char, 2 * seq_size);
  seq_mRNA = generating_mRNA(seq_bin, 0, 2*seq_size);
  //assert_string_equal(expected_char, seq_mRNA);

  expected_char = "AUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG";
  seq_char = "ATCGATCGATCGATCGATCGATCGATCGATCG";
  seq_size = 32;
  char* seq_test = NULL;
  seq_test = calloc(seq_size, sizeof(char));
  char* seq_new = NULL;
  int* ptr;

  for (long int i = 1; i < seq_size; i++) {
    seq_bin = convert_to_binary(seq_char, i);
    seq_new = generating_mRNA(seq_bin, 0, 2 * i);
    ptr = (int*)realloc(seq_test, sizeof(char) * i);
    // Check the realloc worked.
    assert_ptr_not_equal(NULL, ptr);
    // Copy the sequence in the test sequence, to keep only the sequence needed for the test
    memcpy(seq_test, expected_char, i);
    assert_string_equal(seq_test, seq_new);
  }

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_mRNA(NULL, 0, 0));
}


static void test_detecting_genes(void ** state){
  gene_map_t *gene_map = NULL;
  gene_map = malloc(sizeof(*gene_map));

  // Test if the algorithm is OK in a basic case: xxxAUGxxxxUAAxxx
  //                                              AGC AUG AGGGCC UAA CGU
  // The algorithm should detect one gene from the start codon to the stop codon
  
  long int* seq_bin = convert_to_binary("AGCATGAGGGCCTAACGT", 21);
  // printf("sqe_bin : %s \n", seq_bin);
  detecting_genes(seq_bin, 21*2, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(29, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  //                                                         AGC AUG AGGGCC AUG CGAACG UAA CGU
  // The algorithm should detect one gene from the 2nd start codon to the stop codon
  seq_bin = convert_to_binary("AGCATGAGGGCCATGCGAACGTAACGT", 27);
  detecting_genes(seq_bin, 27*2, gene_map);
  assert_int_equal(24, gene_map->gene_start[0]);
  assert_int_equal(47, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  //                                                         AGC AUG CACGCG UAA GCACTG UAA CGU
  // The algorithm should detect one gene from the start codon to the first stop codon
  seq_bin = convert_to_binary("AGCATGCACGCGTAAGCACTGTAACGT", 27);
  detecting_genes(seq_bin, 27*2, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(29, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  //                                                                     CGCCGCGCCGCGGGCG
  // The algorithm should not detect any genes
  gene_map->gene_start[0] = 0;
  gene_map->gene_end[0] = 0;
  seq_bin = convert_to_binary("CGCCGCGCCGCGGGCG", 16);
  detecting_genes(seq_bin, 16*2, gene_map);
  assert_int_equal(0, gene_map->gene_start[0]);
  assert_int_equal(0, gene_map->gene_end[0]);
  assert_int_equal(0, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  //                                                      AGC AUG GCGCAC UAG CGCCCG AUG CUGGGG UAA CGU
  // The algorithm should detect two genes
  seq_bin = convert_to_binary("AGCATGGCGCACTAGCGCCCGATGCTGGGGTAACGT", 36);
  detecting_genes(seq_bin, 36*2, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(29, gene_map->gene_end[0]);
  assert_int_equal(42, gene_map->gene_start[1]);
  assert_int_equal(65, gene_map->gene_end[1]);
  assert_int_equal(2, gene_map->genes_counter);

  free(gene_map->gene_start);
  free(gene_map->gene_end);
  free(gene_map);
}

static void test_generating_aa_chain(void ** state){
  // Test if the algorithm is OK

  // AA -> SEQ -> BINARY -> INVERT BINARY -> CALL FUNCTION
  // K      K      N      N      R      R      S      S      T      T     
  // AAA    AAG    AAC    AAT    AGA    AGG    AGC    AGT    ACA    ACG         
  // 000000 000001 000010 000011 000100 000101 000110 000111 001000 001001   
  assert_string_equal("KKNNRRSSTT", generating_amino_acid_chain((long int []) {0b100100000100111000011000101000001000110000010000100000000000}, 0, 60));
  // T      T      I      M      I      I      E      E      D      D
  // ACC    ACT    ATA    ATG    ATC    ATT    GAA    GAG    GAC    GAT
  // 001010 001011 001100 001101 001110 001111 010000 010001 010010 010011
  assert_string_equal("TTIMIIEEDD", generating_amino_acid_chain((long int []) { 0b110010010010100010000010111100011100101100001100110100010100 }, 0, 60));
  // G      G      A      A      A      A      V      V      V      V 
  // GGC    GGT    GCA    GCG    GCC    GCT    GTA    GTG    GTC    GTT         
  // 010110 010111 011000 011001 011010 011011 011100 011101 011110 011111   
  assert_string_equal("GGAAAAVVVV", generating_amino_acid_chain((long int []) { 0b111110011110101110001110110110010110100110000110111010011010 }, 0, 60));
  // Q      Q      H      H      R      R      R      R      P      P
  // CAA    CAG    CAC    CAT    CGA    CGG    CGC    CGT    CCA    CCG
  // 100000 100001 100010 100011 100100 100101 100110 100111 101000 101001
  assert_string_equal("QQHHRRRRPP", generating_amino_acid_chain((long int []) { 0b100101000101111001011001101001001001110001010001100001000001 }, 0, 60));
  // L      L      L      L      O      O      Y      Y      O      W     
  // CTA    CTG    CTC    CTT    TAA    TAG    TAC    TAT    TGA    TGG    
  // 101100 101101 101110 101111 110000 110001 110010 110011 110100 110101
  assert_string_equal("LLLLOOYYOW", generating_amino_acid_chain((long int []) { 0b101011001011110011010011100011000011111101011101101101001101 }, 0, 60));
  // C      C      S      S      S      S      L      L      F      F
  // TGC    TGT    TCA    TCG    TCC    TCT    TTA    TTG    TTC    TTT
  // 110110 110111 111000 111001 111010 111011 111100 111101 111110 111111
  assert_string_equal("CCSSSSLLFF", generating_amino_acid_chain((long int []) { 0b111111011111101111001111110111010111100111000111111011011011 }, 0, 60));
  // G      G      P      P
  // GGA    GGG    CCC    CCT
  // 010100 010101 101010 101011
  assert_string_equal("GGPP", generating_amino_acid_chain((long int []) { 0b110101010101101010001010 }, 0, 24));


  //  --- Test all the amino acid
  // (alphabetic order of the above sequences.)
  long int* seq_bin = convert_to_binary("AAAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTGAAGAGGACGATGGCGGTGCAGCGGCCGCTGTAGTGGTCGTTCAACAGCACCATCGACGGCGCCGTCCACCGCTACTGCTCCTTTAATAGTACTATTGATGGTGCTGTTCATCGTCCTCTTTATTGTTCTTTGGAGGGCCCCCT", 384);
  char* aa_chain = NULL;
  aa_chain = generating_amino_acid_chain(seq_bin, 0, 384);
  //assert_string_equal("KKNNRRSSTTTTIMIIEEDDGGAAAAVVVVQQHHRRRRPPLLLLOOYYOWCCSSSSLLFFGGPP", aa_chain);

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_amino_acid_chain(NULL, 0, 0));
  // --- invalid value in gene_seq
  // assert_ptr_not_equal("MRGO",
  //                      generating_amino_acid_chain((long int []){1821290092, 18263}, 24));
}


static void test_detecting_mutations(void ** state){
  mutation_map M;
  unsigned short nb_mutations = 6;
  M.size = calloc(nb_mutations, sizeof(unsigned long));
  M.start_mut = calloc(nb_mutations, sizeof(unsigned long));
  M.end_mut = calloc(nb_mutations, sizeof(unsigned long));

  //A : 00
  //T : 11
  //C : 10
  //G : 01

  //GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
  long int* seq_bin = convert_to_binary("GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG", 30);
  //Test if sequence 10 to 23 is a mutation zone and no other mutation zone
  M.size[1]=0;
  M.start_mut[1]=0;
  M.end_mut[1]=0;
  detecting_mutations(seq_bin, 0, 60, M);
  // detecting_mutations((long int []){261725162, 97523700}, 0, 60, M);
  //First mutation is updated with right values
  assert_int_equal(13,M.size[0]);
  assert_int_equal(10,M.start_mut[0]);
  assert_int_equal(24,M.end_mut[0]);
  //No other mutations, should not be updated
  assert_int_equal(0,M.size[1]);
  assert_int_equal(0,M.start_mut[1]);
  assert_int_equal(0,M.end_mut[1]);

  //GTTTTGCAAACGTTAAAGGTTTGAAAGGTG = {261102590, 97523700}
  seq_bin = convert_to_binary("GTTTTGCAAACGTTAAAGGTTTGAAAGGTG", 30);
  //Test if no mutation in this sequence
  M.size[0]=0;
  M.start_mut[0]=0;
  M.end_mut[0]=0;
  detecting_mutations((long int []){261102590, 97523700}, 0, 60, M);
  //No possible mutation zones detected, should not be updated
  assert_int_equal(0,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(0,M.end_mut[0]);

  //GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
  seq_bin = convert_to_binary("GGGCCGTTCCGCCCATAGGCCCGGCTAAGA", 30);
  //Test with 3 mutation zones in this sequence
  M.size[3]=0;
  M.start_mut[3]=0;
  M.end_mut[3]=0;
  detecting_mutations(seq_bin, 0, 60, M);
  //First mutation is updated with right values
  assert_int_equal(11,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(12,M.end_mut[0]);
  //Second mutation is updated with right values
  assert_int_equal(11,M.size[1]);
  assert_int_equal(16,M.start_mut[1]);
  assert_int_equal(28,M.end_mut[1]);
  //Third mutation is updated with right values
  assert_int_equal(15,M.size[2]);
  assert_int_equal(34,M.start_mut[2]);
  assert_int_equal(50,M.end_mut[2]);
  //No other mutations, should not be updated
  assert_int_equal(0,M.size[3]);
  assert_int_equal(0,M.start_mut[3]);
  assert_int_equal(0,M.end_mut[3]);


  // Test three mutations in sequence
  seq_bin = convert_to_binary("GGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTGCGCGCGGCGCGCGGCGCGCGCGCGCGGCGCGCGGCGCGCGGGTTAAAGGTG", 147);
  detecting_mutations(seq_bin, 0, 147*2, M);

  for(int i = 0; i < 3; i++){
    assert_int_equal(85, M.size[i]);
    assert_int_equal(10 + 90*i, M.start_mut[i]);
    assert_int_equal(95 + 90 * i+1, M.end_mut[i]);
  }

  free(M.size);
  free(M.start_mut);
  free(M.end_mut);
}

static void test_calculating_matching_score(void ** state){
  // Test if the algorithm is OK
  // --- With same size
  //  GACCCGAC = 0100101010010010 = {18770}
  //  GGCCAGGC = 0101101000010110 = {26714}
  assert_float_equal(81.250000,
                    calculating_matching_score(
                      (long int []){18770}, 0, 16,
                      (long int []){26714}, 0, 16),
                    0);
  // --- With different size
  //  GACCCGAC = 0100101010010010 = {18770}
  //  TTTCAGGCTC = 11111110000101101110 = {485503}
  assert_float_equal(25.000000,
                    calculating_matching_score(
                      (long int []){18770}, 0, 16,
                      (long int []){485503}, 0, 20),
                    0);
  //  TTTCAGGCTT = 11111110000101101111 = {1009791}
  //  GACCTTCGA = 010010101111100100 = {40786}
  assert_float_equal(45.000000,
                    calculating_matching_score(
                      (long int []){1009791}, 0, 20,
                      (long int []){40786}, 0, 18),
                    0);
  assert_float_equal(35.000000,
                    calculating_matching_score(
                      (long int []){1009791}, 0, 20,
                      (long int []){40786}, 0, 16),
                    0);

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_float_equal(-1.0, calculating_matching_score(NULL, 0, 0, NULL, 0, 0), 0);

}

static void test_get_piece_binary_array(){
  // Test if the algorithm is OK

  long int intsize = int_SIZE + 1;

  long int* arr = (long int []){ 0b1111111111111111111111111111111111111111111111111111111111111111, 0b101000111000101, 0b101000111000101, 0b101000111000101, 0b101000111000101, 31, 31 };
  long int pow = 1;
  for(long int i = 0; i < 64; i ++){
    pow *=2;
    // printf("i : %d => %ld _ %ld\n", i, pow - 1, get_piece_binary_array(arr, 0, i+1)[0]);
    assert_int_equal(pow-1, get_piece_binary_array(arr, 0, i+1)[0]);
  }
  long int* res = get_piece_binary_array(arr, intsize + 1, intsize);
  assert_int_equal(arr[1], res[0]);
  res = get_piece_binary_array(arr, intsize + 1, intsize*5);
  for(long int i =0; i < 5; i++)
    assert_int_equal(arr[i+1], res[i]);

  res = get_piece_binary_array(arr, intsize + 1, 15);
  assert_int_equal(arr[2], res[0]);
  res = get_piece_binary_array(arr, 2*(intsize+1), 15);
  assert_int_equal(arr[3], res[0]);
  res = get_piece_binary_array(arr, 2*(intsize+1), 13);
  assert_int_equal(arr[3]&0b1111111111111, res[0]);

  res = get_piece_binary_array(arr, 2 * intsize, intsize);

  long int size = 31;
  long int pos;

  arr = calloc(size + 1, sizeof(arr));
  for (long int i = 0; i <= size; i++) {
    arr[i] = size - i;
  }

  for (long int it = 0; it <= size; it++) {
    pos = size - it;
    // retourner le Xième bit
    assert_int_equal(pos % 2, get_piece_binary_array(arr, it*(intsize+1), 1)[0]);
    assert_int_equal(pos / 2 % 2, get_piece_binary_array(arr, it*(intsize+1) + 1, 1)[0]);
    assert_int_equal(pos / 4 % 2, get_piece_binary_array(arr, it*(intsize+1) + 2, 1)[0]);
    assert_int_equal(pos / 8 % 2, get_piece_binary_array(arr, it*(intsize+1) + 3, 1)[0]);
    assert_int_equal(pos / 16 % 2, get_piece_binary_array(arr, it*(intsize+1) + 4, 1)[0]);
    // Retourner tout le nombre
    assert_int_equal(pos, get_piece_binary_array(arr, it * (intsize + 1), intsize)[0]);
    // retourner plusieurs bits
    assert_int_equal(pos / 2 % 2 * 2 + pos % 2, get_piece_binary_array(arr, it * (intsize + 1), 2)[0]);
    assert_int_equal(pos / 4 % 2 * 4 + pos / 2 % 2 * 2 + pos % 2, get_piece_binary_array(arr, it * (intsize + 1), 3)[0]);
    assert_int_equal(pos / 4 % 2 * 2 + pos / 2 % 2, get_piece_binary_array(arr, it * (intsize + 1) + 1, 2)[0]);
    assert_int_equal(pos / 16 % 2 * 4 + pos / 8 % 2 * 2 + pos / 4 % 2, get_piece_binary_array(arr, it * (intsize + 1) + 2, 3)[0]);
  }
  free(arr);
}

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
    // BINARIES ARRAYS FUNCTIONS
    cmocka_unit_test(test_get_binary_value),
    // cmocka_unit_test(test_change_binary_value),
    // cmocka_unit_test(test_set_binary_array),
    // cmocka_unit_test(test_xor_binary_array),
    // cmocka_unit_test(test_popcount_binary_array),
    // cmocka_unit_test(test_get_piece_binary_array),
    // // DNA & GENES FUNCTIONS
    // cmocka_unit_test(test_convert_to_binary),
    // cmocka_unit_test(test_binary_to_dna),
    // cmocka_unit_test(test_generating_mRNA),
    // cmocka_unit_test(test_detecting_genes),
    // cmocka_unit_test(test_generating_aa_chain),
    // cmocka_unit_test(test_detecting_mutations),
    // cmocka_unit_test(test_calculating_matching_score),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}
