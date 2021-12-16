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

static void test_calculating_matching_score(void ** state){
  // Test if the algorithm is OK
  // --- With same size
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

int main(void) {
  int result = 0;
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_generating_mRNA),
      cmocka_unit_test(test_detecting_genes),
      cmocka_unit_test(test_generating_aa_chain),
      cmocka_unit_test(test_calculating_matching_score),
      cmocka_unit_test(test_detecting_mutations),
  };
  result |= cmocka_run_group_tests_name("gene", tests, NULL, NULL);

  return result;
}