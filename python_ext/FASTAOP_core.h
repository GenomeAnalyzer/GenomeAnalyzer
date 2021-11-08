#ifndef FASTAOP_CORE_H
#define FASTAOP_CORE_H

int hamming(int b1, int b2);
char* RNA_generate(int seq_len, char* seq);
int genes_detect(int var);
int amino_acid_chain_generate(int var);
int mutations_detect(int var);
int matching_score_calculate(int a, int b);


#endif