#include <stdint.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "gene_bin.h"
#include "gene_bin.c"


int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
    
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}