#include <stdint.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "gene_bin.h"
#include "gene_bin.c"

static void test_get_binary_value(void ** state){
  // Test if the algorithm is OK
      // 1 = 0000000000000000000000000000001
  assert_int_equal(1, get_binary_value((long int []){1}, 0));
      //9350764 = 001101100111010101110001
  assert_int_equal(0, get_binary_value((long int []){1}, 17));
  // Pour chaque binaire de 00000 à 11111, vérifier chaque bit.
  for(int i =0; i<32;i++){
    assert_int_equal(i%2, get_binary_value((long int  []){i}, 0));
    assert_int_equal(i/2%2, get_binary_value((long int  []){i}, 1));
    assert_int_equal(i/4%2, get_binary_value((long int  []){i}, 2));
    assert_int_equal(i/8%2, get_binary_value((long int  []){i}, 3));
    assert_int_equal(i/16%2, get_binary_value((long int  []){i}, 4));
  }
}

static void test_change_binary_value(void ** state){
  // Test if the algorithm is OK

  // 85 = 1010101
  long int  seq_bin = 85;
  int bin;
  // Invert values of seq_bin
  for (int i = 0;i < 7;i++) change_binary_value(&seq_bin, i, i%2);
  // 42 = 0101010
  assert_int_equal(seq_bin, 42);
  // Invert again values of seq_bin.
  for (int i = 0;i < 7;i++) change_binary_value(&seq_bin, i, (i+1)%2);
  assert_int_equal(seq_bin, 85);
}

static void test_set_binary_array(void ** state){
  // Test if the algorithm is OK
  char *seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
  unsigned seq_size = 45;

  long int  *seq_bin = NULL;
  seq_bin = set_binary_array(seq_char, seq_size);

  long int  seq_sol[3] = {2101911378, 172292753, -265825143};

  for (int i = 0; i < 3; ++i)
    assert_int_equal(seq_bin[i], seq_sol[i]);

  free(seq_bin);
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

  // long int  *xor = NULL;
  xor = xor_binary_array(seq_bin, 2 * seq_size, seq_bin2, 2 * seq_size2);
  // xor_size = max size
  int xor_size = seq_size >= seq_size2 ? seq_size : seq_size2;

  long int  xor_sol[3] = {2101911378, 364800657, 2934};

  for (int i = 0; i < 3; ++i)
    assert_int_equal(xor[i], xor_sol[i]);

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

static void test_convert_to_binary(void ** state){
  // Test aa to binary conversions
  // --- Test all valid letters
  long int * res = convert_to_binary("ATCGN", 10);
  assert_int_equal(res[0], 156);


  // Test whether the function correctly detects errors:
  // --- Unknown letter in sequence
  long int * res2 = convert_to_binary("AK", 10);
  assert_ptr_equal(NULL, res2[2]);
}

static void test_binary_to_dna(void ** state){
  // Test binary to aa conversions

  // --- Test all conversion
  assert_string_equal("ATCG", binary_to_dna((long int []){156}, 8));

  // Test whether the function correctly detects errors:
  // --- Wrong size bin_dna_seq
  assert_ptr_equal(NULL, binary_to_dna((long int []){0}, 3));
}

static void test_generating_mRNA(void ** state){
  // Test if the algorithm is OK
      //9350764 = 001101100111010101110001
  assert_string_equal("AUGCGUGGGUAG",
                      generating_mRNA((long int []){9350764}, 24));
      //913666358 = 00110110011101010111000100110110, 30065 = 0111010101110001
  assert_string_equal("AUGCGUGGGUAGAUGCGUGGGUAG",
                      generating_mRNA((long int []){1821290092, 18263}, 48));

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_mRNA(NULL, 0));
  // --- invalid value in gene_seq
      //11957616 = 101101100111010101110000
  assert_ptr_not_equal("AUGCGUGGGUAG",
                       generating_mRNA((long int []){1821290092, 18263}, 48));
}


static void test_detecting_genes(void ** state){
  gene_map_t *gene_map = NULL;
  gene_map = malloc(sizeof(*gene_map));

  // Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the stop codon
  detecting_genes((long int []){963808024, 42}, 40, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(28, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  // The algorithm should detect one gene from the 2nd start codon to the stop codon
  detecting_genes((long int []){732875499, -2036213923}, 60, gene_map);
  assert_int_equal(30, gene_map->gene_start[0]);
  assert_int_equal(48, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the first stop codon
  detecting_genes((long int []){250327787, -2022340747}, 60, gene_map);
  assert_int_equal(10, gene_map->gene_start[0]);
  assert_int_equal(26, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  // The algorithm should not detect any genes
  gene_map->gene_start[0] = 0;
  gene_map->gene_end[0] = 0;
  detecting_genes((long int []){22369621}, 26, gene_map);
  assert_int_equal(0, gene_map->gene_start[0]);
  assert_int_equal(0, gene_map->gene_end[0]);
  assert_int_equal(0, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  // The algorithm should detect two genes
  detecting_genes((long int []) {-469763265, -1612578969, -268435456}, 66, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(24, gene_map->gene_end[0]);
  assert_int_equal(36, gene_map->gene_start[1]);
  assert_int_equal(48, gene_map->gene_end[1]);
  assert_int_equal(2, gene_map->genes_counter);

  free(gene_map->gene_start);
  free(gene_map->gene_end);
  free(gene_map);
}

static void test_generating_aa_chain(void ** state){
  // Test if the algorithm is OK
  //  --- Test all the amino acid
  assert_string_equal("KNKNTTTTRSRSIIIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYSSSSCWCLFLFMOOO",
                      generating_amino_acid_chain((long int []){79823872, -2096862942,
                            -1577991368, 547545866, -1792699787, -1126245655, 1210084514,
                            -752012202, 1001024414, -106443080, -1380064261, -1612777443,
                            189184},
                        384));

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_amino_acid_chain(NULL, 0));
  // --- invalid value in gene_seq
  assert_ptr_not_equal("MRGO",
                       generating_amino_acid_chain((long int []){1821290092, 18263}, 24));
}


static void test_detecting_mutations(void ** state){
  mutation_map M;
  M.size = malloc(5 * sizeof(unsigned long));
  M.start_mut = malloc(5 * sizeof(unsigned long));
  M.end_mut = malloc(5 * sizeof(unsigned long));

  //A : 00
  //T : 11
  //C : 10
  //G : 01

  //GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
  //Test if sequence 10 to 23 is a mutation zone and no other mutation zone
  M.size[1]=0;
  M.start_mut[1]=0;
  M.end_mut[1]=0;
  detecting_mutations((long int []){261725162, 97523700},60, M);
  //First mutation is updated with right values
  assert_int_equal(13,M.size[0]);
  assert_int_equal(10,M.start_mut[0]);
  assert_int_equal(23,M.end_mut[0]);
  //No other mutations, should not be updated
  assert_int_equal(0,M.size[1]);
  assert_int_equal(0,M.start_mut[1]);
  assert_int_equal(0,M.end_mut[1]);

  //GTTTTGCAAACGTTAAAGGTTTGAAAGGTG = {261102590, 97523700}
  //Test if no mutation in this sequence
  M.size[0]=0;
  M.start_mut[0]=0;
  M.end_mut[0]=0;
  detecting_mutations((long int []){261102590, 97523700},60, M);
  //No possible mutation zones detected, should not be updated
  assert_int_equal(0,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(0,M.end_mut[0]);

  //GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
  //Test with 3 mutation zones in this sequence
  M.size[3]=0;
  M.start_mut[3]=0;
  M.end_mut[3]=0;
  detecting_mutations((long int []){-983172758, 17224372},60, M);
  //First mutation is updated with right values
  assert_int_equal(11,M.size[0]);
  assert_int_equal(0,M.start_mut[0]);
  assert_int_equal(11,M.end_mut[0]);
  //First mutation is updated with right values
  assert_int_equal(11,M.size[1]);
  assert_int_equal(16,M.start_mut[1]);
  assert_int_equal(27,M.end_mut[1]);
  //First mutation is updated with right values
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

static void test_calculating_matching_score(void ** state){
  // Test if the algorithm is OK
  // --- With same size
  //  GACCCGAC = 0100101010010010 = {18770}
  //  GGCCAGGC = 0101101000010110 = {26714}
  assert_float_equal(81.250000,
                    calculating_matching_score(
                      (long int []){18770}, 16,
                      (long int []){26714}, 16),
                    0);
  // --- With different size
  //  GACCCGAC = 0100101010010010 = {18770}
  //  TTTCAGGCTC = 11111110000101101110 = {485503}
  assert_float_equal(25.000000,
                    calculating_matching_score(
                      (long int []){18770}, 16,
                      (long int []){485503}, 20),
                    0);
  //  TTTCAGGCTT = 11111110000101101111 = {1009791}
  //  GACCTTCGA = 010010101111100100 = {40786}
  assert_float_equal(45.000000,
                    calculating_matching_score(
                      (long int []){1009791}, 20,
                      (long int []){40786}, 18),
                    0);
  assert_float_equal(35.000000,
                    calculating_matching_score(
                      (long int []){1009791}, 20,
                      (long int []){40786}, 16),
                    0);

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_float_equal(-1.0, calculating_matching_score(NULL, 0, NULL, 0), 0);

}

static void test_get_piece_binary_array(){
  // Test if the algorithm is OK
  long int intsize = int_SIZE + 1;

  long int *arr = (long int []){0b11111111111111111111111111111111, 0b101000111000101, 0b101000111000101, 0b101000111000101, 0b101000111000101, 31, 31};
  long int* res = get_piece_binary_array(arr, 7, 2 * intsize + 1, 2 * intsize + 2);
  assert_int_equal(2, res[0]);
  assert_int_equal(2147494114, res[1]);
  assert_int_equal(2147494114, res[2]);

  res = get_piece_binary_array(arr, 7, 2 * intsize, 2 * intsize + 2);
  assert_int_equal(1, res[0]);
  assert_int_equal(20933, res[1]);
  assert_int_equal(20933, res[2]);

  res = get_piece_binary_array(arr, 7, 2, 3 * intsize + 5);
  assert_int_equal(17, res[0]);
  assert_int_equal(1073747057, res[1]);
  assert_int_equal(1073741831, res[2]);
  assert_int_equal(3221225479, res[3]);

  // Tests récupérations premiers entiers
  assert_int_equal(15, get_piece_binary_array(arr, 7, 0, 4)[0]);
  assert_int_equal(31, get_piece_binary_array(arr, 7, 0, intsize)[0]);

  assert_int_equal(0b101000111000101, get_piece_binary_array(arr, 7, 2 * intsize, intsize)[0]);

  res = get_piece_binary_array(arr, 7, 2 * intsize, intsize);
  assert_int_equal(0b101000111000101, res[0]);
  // assert_int_equal(31, res[1]);

  assert_int_equal(0b11111111111111111111111111111111, get_piece_binary_array(arr, 7, 6 * intsize, intsize)[0]);

  // SURCHARGE DES TESTS AVEC LES TESTS DE mask_binary_array

  // long int size = 511;
  // long int pos;

  // arr = calloc(size + 1, sizeof(arr));
  // for (long int i = 0; i <= size; i++) {
  //   arr[i] = size - i;
  // }

  // for (long int it = 0; it <= size; it++) {
  //   pos = size - it;
  //   // retourner le Xième bit
  //   assert_int_equal(it % 2, get_piece_binary_array(arr, size+1, 0, 1)[0]);
  //   assert_int_equal(it / 2 % 2, get_piece_binary_array(arr, size+1, 1, 1)[0]);
  //   assert_int_equal(it / 4 % 2, get_piece_binary_array(arr, size+1, 2, 1)[0]);
  //   assert_int_equal(it / 8 % 2, get_piece_binary_array(arr, size+1, 3, 1)[0]);
  //   assert_int_equal(it / 16 % 2, get_piece_binary_array(arr, size+1, 4, 1)[0]);
  //   // Retourner tout le nombre
  //   assert_int_equal(it, get_piece_binary_array(arr, size+1, 0, intsize)[0]);
  //   // retourner plusieurs bits
  //   assert_int_equal(it / 2 % 2 * 2 + it % 2, get_piece_binary_array(arr, size+1, 0, 2)[0]);
  //   assert_int_equal(it / 4 % 2 * 4 + it / 2 % 2 * 2 + it % 2, get_piece_binary_array(arr, size+1, 0, 3)[0]);
  //   assert_int_equal(it / 4 % 2 * 2 + it / 2 % 2, get_piece_binary_array(arr, size+1, 1, 2)[0]);
  //   assert_int_equal(it / 16 % 2 * 4 + it / 8 % 2 * 2 + it / 4 % 2, get_piece_binary_array(arr, size+1, 2, 3)[0]);
  // }
  // free(arr);
}

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
    // BINARIES ARRAYS FUNCTIONS
    cmocka_unit_test(test_get_binary_value),
    cmocka_unit_test(test_change_binary_value),
    cmocka_unit_test(test_set_binary_array),
    cmocka_unit_test(test_xor_binary_array),
    cmocka_unit_test(test_popcount_binary_array),
    cmocka_unit_test(test_get_piece_binary_array),
    // DNA & GENES FUNCTIONS
    cmocka_unit_test(test_convert_to_binary),
    cmocka_unit_test(test_binary_to_dna),
    cmocka_unit_test(test_generating_mRNA),
    cmocka_unit_test(test_detecting_genes),
    cmocka_unit_test(test_generating_aa_chain),
    cmocka_unit_test(test_detecting_mutations),
    cmocka_unit_test(test_calculating_matching_score),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}