#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include "rdtsc.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

void print_sim_mat_color(int* F, char A [], char B [], int m, int n, int k, int l) {
    int couleur = 81;
    printf("%2s %2s ", "*", "-");
    for (int j = 0; j < m; j++)
        printf("%2c ", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%2c ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 2d \033[0m", F[i * m + j]);
            else if (i >= 1 && j >= 1) {
                printf("\033[29m\033[1m\033[48;5;%dm% 2d \033[0m", couleur, F[i * m + j]);
                couleur--;
            }
            else
                printf("% 2d ", F[i * m + j]);
        }
        printf("\n");
    }
}

void print_sim_mat(int* F, char A [], char B [], int m, int n, int k, int l) {
    printf(" \033[48;5;8m%2s  %3s  \033[0m", "*", "-");
    for (int j = 0; j < m - 1; j++)
        printf("\033[48;5;8m%3c  \033[0m", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf(" \033[48;5;8m%2c \033[0m ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 3d  \033[0m", F[i * n + j]);
            else
                printf("% 3d  ", F[i * n + j]);
        }
        printf("\n");
    }
}

int* calculate_scoring_matrix(char A [], char B [], int match, int mismatch, int gap) {

    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < m; i++)
        F[i] = gap * i;
    for (int j = 0; j < n; j++)
        F[j * m] = gap * j;

    // Calculate scoring matrix
    for (int i = 1; i < n; i++)
        for (int j = 1; j < m; j++) {
            diag = F[(i - 1) * m + j - 1] + (A[j - 1] == B[i - 1] ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    return F;
}

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




int main() {
    // char A [] = "GCATGCTATA";
    // char B [] = "GATTGCATAG";

    char A [] = "GCATGC";
    char B [] = "GATTA";

    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int match = 1;
    int mismatch = -1;
    int gap = -1;

    unsigned long long before, after;
    double elapsed;

    before = rdtsc();
    int* F = calculate_scoring_matrix(A, B, match, mismatch, gap);
    after = rdtsc();
    print_sim_mat(F, A, B, m, n, n - 1, m - 1);
    // print_sim_mat_color(F, A, B, m, n, n - 1, m - 1);

    // int score = align(F, A, B, match, mismatch, gap, 1);
    int score = F[m * n - 1];
    elapsed = (double)(after - before);
    printf("*Max score : %d\n", score);
    printf("*Cycles count : %.0f\n", elapsed);
    printf("*Time : %.3f sec\n", elapsed / 3.3e9);
    free(F);
}