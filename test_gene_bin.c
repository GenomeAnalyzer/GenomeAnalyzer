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
  assert_int_equal(1, get_binary_value((unsigned int[]){1}, 0));
      //9350764 = 001101100111010101110001
  assert_int_equal(0, get_binary_value((unsigned int[]){1}, 17));
}

static void test_change_binary_value(void ** state){
  // Test if the algorithm is OK
      //9350764 = 001101100111010101110001
  unsigned int seq_bin = 9350764;
  int bin_5 = get_binary_value(&seq_bin, 5);
  change_binary_value(&seq_bin, 5, 0);
  assert_int_equal(0, get_binary_value(&seq_bin, 5));
  assert_int_not_equal(bin_5, get_binary_value(&seq_bin, 5));
}

static void test_set_binary_array(void ** state){
  // Test if the algorithm is OK
  char *seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
  unsigned seq_size = 45;

  unsigned int *seq_bin = NULL;
  seq_bin = set_binary_array(seq_char, seq_size);

  unsigned int seq_sol[3] = {2101911378, 172292753, -265825143};

  for (int i = 0; i < 3; ++i)
    assert_int_equal(seq_bin[i], seq_sol[i]);

  free(seq_bin);
}

static void test_xor_binary_array(void ** state){
  // Test if the algorithm is OK
  char *seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
  unsigned seq_size = 45;
  unsigned int *seq_bin = NULL;
  seq_bin = set_binary_array(seq_char, seq_size);

  char *seq_char2 = "GACCTTTTTTTTTTTTTCTTCGA";
  unsigned seq_size2 = 23;
  unsigned int *seq_bin2 = NULL;
  seq_bin2 = set_binary_array(seq_char2, seq_size2);

  unsigned int *xor = NULL;
  xor = xor_binary_array(seq_bin, 2 * seq_size, seq_bin2, 2 * seq_size2);
  // xor_size = max size
  int xor_size = seq_size >= seq_size2 ? seq_size : seq_size2;

  unsigned int xor_sol[3] = {2101911378, 364800657, 2934};

  for (int i = 0; i < 3; ++i)
    assert_int_equal(xor[i], xor_sol[i]);

  free(seq_bin);
  free(seq_bin2);
  free(xor);
}

static void test_popcount_binary_array(void ** state){
  // Test if the algorithm is OK
      // 1 = 0000000000000000000000000000001
  assert_int_equal(1, popcount_binary_array((unsigned int[]){1}, 31));
      //9350764 = 001101100111010101110001
  assert_int_equal(13, popcount_binary_array((unsigned int[]){9350764}, 24));
}

static void test_convert_to_binary(void ** state){
  // Test aa to binary conversions
  // --- Test all valid letters
  unsigned int* res = convert_to_binary("ATCGN", 10);
  assert_int_equal(res[0], 156);


  // Test whether the function correctly detects errors:
  // --- Unknown letter in sequence
  unsigned int* res2 = convert_to_binary("AK", 10);
  assert_ptr_equal(NULL, res2[2]);
}

static void test_binary_to_dna(void ** state){
  // Test binary to aa conversions

  // --- Test all conversion
  assert_string_equal("ATCG", binary_to_dna((unsigned int[]){156}, 8));

  // Test whether the function correctly detects errors:
  // --- Wrong size bin_dna_seq
  assert_ptr_equal(NULL, binary_to_dna((unsigned int[]){0}, 3));
}

static void test_generating_mRNA(void ** state){
  // Test if the algorithm is OK
      //9350764 = 001101100111010101110001
  assert_string_equal("AUGCGUGGGUAG",
                      generating_mRNA((unsigned int[]){9350764}, 24));
      //913666358 = 00110110011101010111000100110110, 30065 = 0111010101110001
  assert_string_equal("AUGCGUGGGUAGAUGCGUGGGUAG",
                      generating_mRNA((unsigned int[]){1821290092, 18263}, 48));

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_mRNA(NULL, 0));
  // --- invalid value in gene_seq
      //11957616 = 101101100111010101110000
  assert_ptr_not_equal("AUGCGUGGGUAG",
                       generating_mRNA((unsigned int[]){1821290092, 18263}, 48));
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
    // DNA & GENES FUNCTIONS
    cmocka_unit_test(test_convert_to_binary),
    cmocka_unit_test(test_binary_to_dna),
    cmocka_unit_test(test_generating_mRNA),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}