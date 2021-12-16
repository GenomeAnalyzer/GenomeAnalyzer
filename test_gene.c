#include <stdint.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

#include "gene.c"
#include "gene.h"

static void test_generating_mRNA(void ** state){
  // Test if the algorithm is OK
  assert_string_equal("AUGCGUGGGUAG",
                      generating_mRNA((unsigned short[]){0,0,1,1,0,1,1,0,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,1},
                                      24));

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_mRNA(NULL, 0));
  // --- invalid value in gene_seq
  assert_ptr_not_equal("AUGCGUGGGUAG",
                       generating_mRNA((unsigned short[]){0,0,1,44,0,1,1,0,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,1},
                                       24));
}

static void test_detecting_genes(void ** state){
  gene_map_t *gene_map = NULL;
  gene_map = malloc(sizeof(*gene_map));

  // Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the stop codon
  detecting_genes((unsigned int[])
    {0,0,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,0},
                  40, gene_map);
  assert_int_equal(6, gene_map->gene_start[0]);
  assert_int_equal(28, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  // The algorithm should detect one gene from the 2nd start codon to the stop codon
  detecting_genes((unsigned int[])
    {1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1},
                  60, gene_map);
  assert_int_equal(30, gene_map->gene_start[0]);
  assert_int_equal(48, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  // The algorithm should detect one gene from the start codon to the first stop codon
  detecting_genes((unsigned int[])
    {1,1,0,1,0,1,1,1,0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,1,1},
                  60, gene_map);
  assert_int_equal(10, gene_map->gene_start[0]);
  assert_int_equal(26, gene_map->gene_end[0]);
  assert_int_equal(1, gene_map->genes_counter);

  // Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  // The algorithm should not detect any genes
  gene_map->gene_start[0] = 0;
  gene_map->gene_end[0] = 0;
  detecting_genes((unsigned int[])
    {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0},
                  26, gene_map);
  assert_int_equal(0, gene_map->gene_start[0]);
  assert_int_equal(0, gene_map->gene_end[0]);
  assert_int_equal(0, gene_map->genes_counter);

  // Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  // The algorithm should detect two genes
  detecting_genes((unsigned int[])
    {1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1},
                  66, gene_map);
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
                      generating_amino_acid_chain((unsigned short[]){0,0,0,0,0,0,0,0,0,0,1,0,
                         0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,
                         1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,
                         1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,1,1,
                         1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,0,0,1,0,0,1,0,0,1,
                         1,0,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,1,0,1,1,0,1,1,0,
                         1,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,1,0,1,1,0,0,0,
                         0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,
                         0,1,0,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,1,1,
                         0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,
                         1,1,0,1,1,0,1,1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,1,
                         0,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0},
                        384));

  // Test whether the function correctly detects errors:
  // --- NULL error
  assert_ptr_equal(NULL, generating_amino_acid_chain(NULL, 0));
  // --- invalid value in gene_seq
  assert_ptr_not_equal("MRGO",
                       generating_amino_acid_chain((unsigned short[]){0,0,1,44,0,1,1,0,0,1,1,1,
                        0,1,0,1,0,1,1,1,0,0,0,1},
                       24));

}

//Tests for detecting_mutations function
static void test_detecting_mutations(void ** state){
  //The function should return true if we give in input an array with a 0,1 or 1,0 continuous sequence
  //higher than 1/5th of the half size of the array.
 assert_true(detecting_mutations((unsigned int[]){0,1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1},60));
  //The function should return false if we give in input an array with a 0,1 or 1,0 continuous sequence
  //lower than 1/5th of the half size of the array.
  assert_false(detecting_mutations((unsigned int[]){0,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1},60));
}

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_generating_mRNA),
      cmocka_unit_test(test_detecting_genes),
      cmocka_unit_test(test_generating_aa_chain),
      cmocka_unit_test(test_detecting_mutations),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}