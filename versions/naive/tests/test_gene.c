#include <stdint.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>
#include <string.h>

#include "../src/gene.c"
#include "../headers/gene.h"

static void test_generating_mRNA(void ** state){
  // Test if the algorithm is OK

  unsigned short bin_dna_seq[16] = { 0,0,1,1,1,0,0,1 };
  assert_string_equal("AUCG", generating_mRNA(bin_dna_seq, 8));

  unsigned short bin_dna_seq2[16] = { 0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1 };
  assert_string_equal("AUCGAUCG", generating_mRNA(bin_dna_seq2, 16));
  assert_string_equal("AUCGAUC", generating_mRNA(bin_dna_seq2, 14));

  char* seq_char = "AUCGAUCGAUCGAUCG";
  unsigned short seq_size = 16;
  unsigned short bin_dna_seq4[32] = { 0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1 };
  char* seq_test = NULL;
  seq_test = calloc(seq_size, sizeof(char));
  char* seq_new = NULL;
  int* ptr;

  for (long int i = 1; i < seq_size; i++) {
    seq_new = generating_mRNA(bin_dna_seq4, 2 * i);
    ptr = (int*)realloc(seq_test, sizeof(char) * i);
    // Check the realloc worked.
    assert_ptr_not_equal(NULL, ptr);
    // Copy the sequence in the test sequence, to keep only the sequence needed for the test
    memcpy(seq_test, seq_char, i);
    assert_string_equal(seq_test, seq_new);
  }
}

static void test_detecting_genes(void ** state){
  gene_map_t *gene_map = NULL;
  gene_map = malloc(sizeof(*gene_map));

  // Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the stop codon
  detecting_genes((unsigned long[])
    {0,0,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,0},
                  40, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(33, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  // The algorithm should detect one gene from the 2nd start codon to the stop codon
  detecting_genes((unsigned long[])
    {1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1},
                  60, gene_map);
  assert_int_equal(30, gene_map->gene_start[0]);
  assert_int_equal(53, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the first stop codon
  detecting_genes((unsigned long[])
    {1,1,0,1,0,1,1,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,1,1},
                  60, gene_map);
  assert_int_equal(10, gene_map->gene_start[0]);
  assert_int_equal(31, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  // The algorithm should not detect any genes
  gene_map->gene_start[0] = 0;
  gene_map->gene_end[0] = 0;
  detecting_genes((unsigned long[])
    {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0},
                  26, gene_map);
  assert_int_equal(0, gene_map->gene_start[0]);
  assert_int_equal(0, gene_map->gene_end[0]);
  assert_int_equal(0, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  // The algorithm should detect two genes
  detecting_genes((unsigned long[])
    {1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1},
                  66, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(29, gene_map->gene_end[0]);
  assert_int_equal(36, gene_map->gene_start[1]);
  assert_int_equal(53, gene_map->gene_end[1]);
  assert_int_equal(2, gene_map->genes_counter);

  free(gene_map->gene_start);
  free(gene_map->gene_end);
  free(gene_map);
}

static void test_generating_aa_chain(void ** state){
  // Test if the algorithm is OK

  // AA -> SEQ -> BINARY -> CALL FUNCTION
  // K      K      N      N      R      R      S      S      T      T     
  // AAA    AAG    AAC    AAT    AGA    AGG    AGC    AGT    ACA    ACG         
  // 000000 000001 000010 000011 000100 000101 000110 000111 001000 001001
  assert_string_equal("KKNNRRSSTT", generating_amino_acid_chain((unsigned short[]) {0,0,0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1 }, 60));
  // T      T      I      M      I      I      E      E      D      D
  // ACC    ACT    ATA    ATG    ATC    ATT    GAA    GAG    GAC    GAT
  // 001010 001011 001100 001101 001110 001111 010000 010001 010010 010011
  assert_string_equal("TTIMIIEEDD", generating_amino_acid_chain((unsigned short[]) {0,0,1,0,1,0,0,0,1,0,1,1,0,0,
  1,1,0,0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,1 }, 60));
  // G      G      A      A      A      A      V      V      V      V 
  // GGC    GGT    GCA    GCG    GCC    GCT    GTA    GTG    GTC    GTT         
  // 010110 010111 011000 011001 011010 011011 011100 011101 011110 011111
  assert_string_equal("GGAAAAVVVV", generating_amino_acid_chain((unsigned short[]) {0,1,0,1,1,0,0,1,0,1,1,1,0,1,
  1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,1,1 }, 60));
  // Q      Q      H      H      R      R      R      R      P      P
  // CAA    CAG    CAC    CAT    CGA    CGG    CGC    CGT    CCA    CCG
  // 100000 100001 100010 100011 100100 100101 100110 100111 101000 101001
  assert_string_equal("QQHHRRRRPP", generating_amino_acid_chain((unsigned short[]) {1,0,0,0,0,0,1,0,0,0,0,1,1,0,
  0,0,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0,1,0,0,1 }, 60));
  // L      L      L      L      O      O      Y      Y      O      W     
  // CTA    CTG    CTC    CTT    TAA    TAG    TAC    TAT    TGA    TGG    
  // 101100 101101 101110 101111 110000 110001 110010 110011 110100 110101
  assert_string_equal("LLLLOOYYOW", generating_amino_acid_chain((unsigned short[]) {1,0,1,1,0,0,1,0,1,1,0,1,1,0,
  1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,1,0,1 }, 60));
  // C      C      S      S      S      S      L      L      F      F
  // TGC    TGT    TCA    TCG    TCC    TCT    TTA    TTG    TTC    TTT
  // 110110 110111 111000 111001 111010 111011 111100 111101 111110 111111
  assert_string_equal("CCSSSSLLFF", generating_amino_acid_chain((unsigned short[]) {1,1,0,1,1,0,1,1,0,1,1,1,1,1,
  1,0,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1 }, 60));
  // G      G      P      P
  // GGA    GGG    CCC    CCT
  // 010100 010101 101010 101011
  assert_string_equal("GGPP", generating_amino_acid_chain((unsigned short []) { 0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,1 }, 24));


  //  --- Test all the amino acid
  // (alphabetic order of the above sequences.)
  unsigned short* seq_bin = convert_to_binary("AAAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTGAAGAGGACGATGGCGGTGCAGCGGCCGCTGTAGTGGTCGTTCAACAGCACCATCGACGGCGCCGTCCACCGCTACTGCTCCTTTAATAGTACTATTGATGGTGCTGTTCATCGTCCTCTTTATTGTTCTTTGGAGGGCCCCCT", 384);
  char* aa_chain = NULL;
  aa_chain = generating_amino_acid_chain(seq_bin, 384);
  assert_string_equal("KKNNRRSSTTTTIMIIEEDDGGAAAAVVVVQQHHRRRRPPLLLLOOYYOWCCSSSSLLFFGGPP", aa_chain);

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_amino_acid_chain(NULL, 0));
  // --- invalid value in gene_seq
  assert_ptr_not_equal("MRGO",
                       generating_amino_acid_chain((unsigned short[]){0,0,1,44,0,1,1,0,0,1,1,1,
                        0,1,0,1,0,1,1,1,0,0,0,1},
                       24));
}

static void test_convert_to_binary(void ** state){
  // Test aa to binary conversions
  unsigned short* res = convert_to_binary("ATCG", 8);
  unsigned short expected_res = 0b0011100100;

  // --- Test all valid letters
  // A 00 T 11 C 10 G 01
  for(unsigned short i = 0; i < 10; i++)
    assert_int_equal((expected_res >> (9-i)) & 0b1, res[i]);

  expected_res = 0b0011100100111001;
  res = convert_to_binary("ATCGATCG", 16);
  for(unsigned short i = 0; i < 16; i++)
    assert_int_equal((expected_res >> (15-i)) & 0b1, res[i]);

  // Give a different size than the char sequence.
  res = convert_to_binary("ATCGATCG", 8);
  expected_res = 0b00111001;
  for (unsigned short i = 0; i < 8; i++)
    assert_int_equal((expected_res >> (7 - i)) & 0b1, res[i]);

  res = convert_to_binary("ATCGATCG", 7);
  assert_int_equal((expected_res >> 7) & 0b1, res[0]);
  for (unsigned short i = 0; i < 7; i++)
    assert_int_equal((expected_res >> (7 - i)) & 0b1, res[i]);

  // Test whether the function correctly detects errors:
  unsigned short* res2 = convert_to_binary("AK", 10);

  // --- Unknown letter in sequence
  assert_ptr_equal(NULL, res2);
}

static void test_binary_to_aa(void ** state){
  // Test binary to aa conversions

  // --- Test all conversion
  unsigned short bin_dna_seq[16] = { 0,0,1,1,1,0,0,1 };
  assert_string_equal("ATCG", binary_to_aa(bin_dna_seq, 8));

  // Test whether the function correctly detects errors:
  // --- Wrong size bin_dna_seq
  unsigned short bin_dna_seq2[3] = {0,0,1};
  assert_ptr_equal(NULL, binary_to_aa(bin_dna_seq2, 3));

  unsigned short bin_dna_seq3[16] = {0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1};
  assert_string_equal("ATCGATCG", binary_to_aa(bin_dna_seq3, 16));
  assert_string_equal("ATCGATC", binary_to_aa(bin_dna_seq3, 14));

  char* seq_char = "ATCGATCGATCGATCG";
  unsigned short seq_size = 16;
  unsigned short bin_dna_seq4[32] = { 0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1 };
  char* seq_test = NULL;
  seq_test = calloc(seq_size, sizeof(char));
  char* seq_new = NULL;
  int* ptr;

  for (long int i = 1; i < seq_size; i++) {
    seq_new = binary_to_aa(bin_dna_seq4, 2 * i);
    ptr = (int*)realloc(seq_test, sizeof(char) * i);
    // Check the realloc worked.
    assert_ptr_not_equal(NULL, ptr);
    // Copy the sequence in the test sequence, to keep only the sequence needed for the test
    memcpy(seq_test, seq_char, i);
    assert_string_equal(seq_test, seq_new);
  }

}

static void test_calculating_matching_score(void ** state){
  // Test if the algorithm is OK
  //A : 00 T : 11 C : 10 G : 01

  // --- With same size
  // ATCGT = 0011100111
  unsigned short seq_bin[10] = {0,0,1,1,1,0,0,1,1,1}; 
  // Check exact same sequence and size
  assert_float_equal(100, calculating_matching_score(seq_bin, 10, seq_bin, 10), 0);
  // Check same sequence but size null for second parameter.
  assert_float_equal(40, calculating_matching_score(seq_bin, 10, seq_bin, 0), 0);

  assert_float_equal(81.250000,
                    calculating_matching_score(
                      (unsigned short[]){0,1,0,0,1,0,1,0,1,0,0,1,0,0,1,0}, 16,
                      (unsigned short[]){0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0}, 16),
                    0);
  // --- With different size
  assert_float_equal(42.105263,
                    calculating_matching_score(
                      (unsigned short[]){0,1,0,0,1,0,1,0,1,0,0,1,0,0,1,0}, 16,
                      (unsigned short[]){1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1}, 19),
                    0);

  assert_float_equal(47.368420,
                    calculating_matching_score(
                      (unsigned short[]){1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,1,1,1}, 19,
                      (unsigned short[]){0,1,0,0,1,0,1,0,1,1,1,1,1,0,0,1,0}, 17),
                    0);

  // --- Test the matching score for each size
  unsigned short ar1[20] = { 0 };
  unsigned short ar2[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (short i = 0; i < 20; i++) {
    assert_float_equal( (float)(i) * 5,
                        calculating_matching_score(ar1, i + 1, ar2, 20), 0);
    ar1[i] = 1;
    assert_float_equal( (float)(i + 1) * 5,
                        calculating_matching_score(ar1, i + 1, ar2, 20), 0);
  }

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_float_equal(-1.0, calculating_matching_score(NULL, 0, NULL, 0), 0);

}

//Tests for detecting_mutations function
static void test_detecting_mutations(void ** state){
  mutation_map M;
  M.size = malloc(5 * sizeof(unsigned long));
  M.start_mut = malloc(5 * sizeof(unsigned long));
  M.end_mut = malloc(5 * sizeof(unsigned long));

  //A : 00
  //T : 11
  //C : 10
  //G : 01

  //GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG
  //Test if sequence 10 to 23 is a mutation zone and no other mutation zone
  M.size[1]=0;
  M.start_mut[1]=0;
  M.end_mut[1]=0;
  detecting_mutations((unsigned short[]){0,1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1},60, M);
  //First mutation is updated with right values
  assert_int_equal(13,M.size[0]);
  assert_int_equal(10,M.start_mut[0]);
  assert_int_equal(23,M.end_mut[0]);
  //No other mutations, should not be updated
  assert_int_equal(0,M.size[1]);
  assert_int_equal(0,M.start_mut[1]);
  assert_int_equal(0,M.end_mut[1]);

  //GTTTTGCAAACGTTAAAGGTTTGAAAGGTG
  //Test if no mutation in this sequence
  M.size[0]=0;
  M.start_mut[0]=0;
  M.end_mut[0]=0;
  detecting_mutations((unsigned short[]){0,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1},60, M);
  //No possible mutation zones detected, should not be updated
  assert_int_equal(0,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(0,M.end_mut[0]);

  //GGGCCGTTCCGCCCATAGGCCCGGCTAAGA
  //Test with 3 mutation zones in this sequence
  M.size[3]=0;
  M.start_mut[3]=0;
  M.end_mut[3]=0;
  detecting_mutations((unsigned short[]){0,1,0,1,0,1,1,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1,1,0,1,0,1,0,0,0,1,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0},60, M);
  //First mutation is updated with right values
  assert_int_equal(11,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(11,M.end_mut[0]);
  //Second mutation is updated with right values
  assert_int_equal(11,M.size[1]);
  assert_int_equal(16,M.start_mut[1]);
  assert_int_equal(27,M.end_mut[1]);
  //Third mutation is updated with right values
  assert_int_equal(15,M.size[2]);
  assert_int_equal(34,M.start_mut[2]);
  assert_int_equal(49,M.end_mut[2]);
  //No other mutations, should not be updated
  assert_int_equal(0,M.size[3]);
  assert_int_equal(0,M.start_mut[3]);
  assert_int_equal(0,M.end_mut[3]);

  free(M.size);
  free(M.start_mut);
  free(M.end_mut);
  }

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_generating_mRNA),
      cmocka_unit_test(test_detecting_genes),
      cmocka_unit_test(test_generating_aa_chain),
      cmocka_unit_test(test_calculating_matching_score),
      cmocka_unit_test(test_detecting_mutations),
      cmocka_unit_test(test_convert_to_binary),
      cmocka_unit_test(test_binary_to_aa),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}
