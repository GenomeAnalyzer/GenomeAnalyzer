#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "rdtsc.h"

#include "gene_bin.h"
#include "gene_bin.c"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

//////////////// Print Similarity Matrix
/**
 * Prints the similarity matrix
 *
 * in : F : The similarity matrix in COL MAJOR. Its dimensions are `m*n`.
 * in : A : Horizontal char array. Its size is `m-1`
 * in : B : Vertical char array. Its size is `n-1`
 * in : m : Number of columns of `F`. Also the `A` length + 1
 * in : n : Number of lines of `F`. Also the `B` length + 1
 * in : k : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 * in : l : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 */
void print_sim_mat(int* F, char A [], char B [], int m, int n, int k, int l) {
    printf("%3s  %3s  ", "*", "-");
    for (int j = 0; j < m; j++)
        printf("%3c  ", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%3c  ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 3d  \033[0m", F[i * m + j]);
            else
                printf("% 3d  ", F[i * m + j]);
        }
        printf("\n");
    }
}


//////////////// Calculate Similarity Matrix
/**
 * Calculates the similarity matrix of two sequences.
 * Sequences are stored as pairs of two bits, for each nucleotide.
 * One extra space for the similarity matrix is required, its returned dimensions are `(m+1)*(n+1)`
 *
 * in : bin_A : Horizontal binary array. Its size is `m*2`
 * in : bin_B : Vertical binary array. Its size is `n*2`
 * in : m : Number of columns of `F`. Also the `A` length + 1
 * in : n : Number of lines of `F`. Also the `B` length + 1
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * out : F : The similarity matrix in COL MAJOR. <!> Its dimensions are `(m+1)*(n+1)`.
 */
int* bin_calculate_scoring_matrix(long int* bin_A, long int* bin_B, int m, int n, int match, int mismatch, int gap, char A [], char B []) {
    m++;
    n++;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < n; i++)
        F[i] = gap * i;
    for (int j = 0; j < m; j++)
        F[j * n] = gap * j;

    // Relative position of iterator in bin_X
    int i_bin_A, j_bin_B;

    // Two right bits
    unsigned short int aa, bb;

    //    G C A T G C G 
    //<=> G C G T A C G 
    //  0b10011011000110

    // Calculate scoring matrix
    for (int i = 1; i < n; i++) {
        i_bin_A = i / int_SIZE + (i % int_SIZE != 0) - 1;
        aa = (bin_A[i_bin_A] >> ((i - 1) * 2) % int_SIZE) & 0b11;
        for (int j = 1; j < m; j++) {
            j_bin_B = j / int_SIZE + (j % int_SIZE != 0) - 1;
            bb = (bin_B[j_bin_B] >> ((j - 1) * 2) % int_SIZE) & 0b11;
            diag = F[(i - 1) * m + j - 1] + (aa == bb ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    }
    return F;
}

//////////////// String Insert
/**
 * Insert a char at beginning of a char array.
 *
 * in/out : dest   : Destination string, to which we insert the `source` char
 * in     : source : Char to be inerted at beginning of `dest`.
 */
void string_insert(char** dest, char* source) {
    if (*dest == NULL) {
        *dest = calloc(sizeof(char), 2);
        memcpy(*dest, source, 1);
        return;
    }
    char* temp = calloc(sizeof(char), strlen(*dest) + 2);
    memcpy(temp, source, 1);
    strcat(temp, *dest);
    free(*dest);
    *dest = temp;
}

//////////////// Align
/**
 * Backtrace the similarity matrix, and prints the two sequences aligned.
 * Returns the max similarity score of the two sequences.
 *
 * in  : F : The similarity matrix in COL MAJOR. Its dimensions are given by (`A` length + 1)*(`B` length + 1).
 * in  : A : Horizontal char array.
 * in  : B : Vertical char array.
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * in  : print : Enables (1) or Disables (0) alignment output. Returns max similarity score either ways.
 * out : Max similarity score of the `A` and `B` sequences.
 */
int align(int* F, char A [], char B [], int match, int mismatch, int gap, int print) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    if (print) {
        int i = m - 1;
        int j = n - 1;

        char* A_aligned = NULL;
        char* B_aligned = NULL;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && F[i * m + j] == F[(i - 1) * m + j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch)) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, B + (j - 1));
                i--;
                j--;
            }
            else if (i > 0 && F[i * m + j] == F[(i - 1) * m + j] + gap) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, "-");
                i--;
            }
            else {
                // else if (j > 0 && F[i * m + j] == F[i * m + j - 1] + gap) {
                string_insert(&A_aligned, "-");
                string_insert(&B_aligned, B + (j - 1));
                j--;
            }
        }
        printf("%s\n", A_aligned);
        printf("%s\n", B_aligned);
    }
    return F[m * n - 1];
}



// Basic main test, using bin's functions.
int main() {
    // 0b10011011000110
    char A [] = "GCATGCG";
    // 0b00010011110010
    char B [] = "GATTACAAA";

    long int* bin_A = set_binary_array(A, strlen(A));
    long int* bin_B = set_binary_array(B, strlen(B));

    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int match = 1;
    int mismatch = -1;
    int gap = -1;

    unsigned long long before, after;
    double elapsed;

    // int* F = calculate_scoring_matrix(A, B, match, mismatch, gap);
    before = rdtsc();
    int* F = bin_calculate_scoring_matrix(bin_A, bin_B, strlen(A), strlen(B), match, mismatch, gap, A, B);
    after = rdtsc();
    // print_sim_mat(F, A, B, m, n, n - 1, m - 1);

    // int score = align(F, A, B, match, mismatch, gap, 0);
    int score = F[m * n - 1];

    elapsed = (double)(after - before);
    printf("*Max score : %d\n", score);
    printf("*Cycles count : %.0f\n", elapsed);
    printf("*Time : %.3f sec\n", elapsed / 3.3e9);
    free(F);
}