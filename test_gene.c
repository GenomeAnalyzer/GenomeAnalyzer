#include <stdint.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "gene.c"

//Tests for detecting_mutations function
static void test_detecting_mutations(void ** state){

  //The function should return true if we give in input an array with a 0,1 or 1,0 continuous sequence
  //higher than 1/5th of the half size of the array.
  assert_true(detecting_mutations((int[]){0,1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,1,0,0
                                    ,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1}, 60));
  //The function should return false if we give in input an array with a 0,1 or 1,0 continuous sequence
  //lower than 1/5th of the half size of the array.
  assert_false(detecting_mutations((int[]){0,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0
                                    ,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1}, 60));
}

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_detecting_mutations),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}