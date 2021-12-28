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

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_get_binary_value),
    cmocka_unit_test(test_change_binary_value),
    cmocka_unit_test(test_set_binary_array),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}