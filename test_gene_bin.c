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

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_get_binary_value),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}