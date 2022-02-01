#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x, y) ((x) > (y) ? (x) : (y))

void print_mat(int* F, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d\t", F[i * n + j]);
        }
        printf("\n");
    }
}

int* calculate_scoring_matrix(char A [], char B [], int match, int mismatch, int gap) {

    int m = strlen(A) + 1;
    int n = strlen(B) + 1;
    // printf("m:%d _ n:%d\n", m, n);

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int test = 0;
    int left = 0;
    int up = 0;

    // Initialize F
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            F[i * n + j] = 0;

    // Setting F borders
    for (int i = 0; i < n; i++)
        F[i] = gap * i;
    for (int j = 0; j < m; j++)
        F[j * n] = gap * j;

    // Calculate scoring matrix
    for (int i = 1; i < m; i++)
        for (int j = 1; j < n; j++) {
            test = F[(i - 1) * n + j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch);
            left = F[(i - 1) * n + j] + gap;
            up = F[i * n + j - 1] + gap;
            int var = MAX(MAX(test, left), up);
            F[i * n + j] = var;
        }
    return F;
}

void string_insert(char* dest, char* source) {
    // Insert first source char into dest.
    dest = realloc(dest, sizeof(dest) + sizeof(char));
    strncpy(dest + 1, dest, strlen(dest));
    strncpy(dest, source, 1);
}

void align(int* F, char A [], char B [], int match, int mismatch, int gap) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    printf("%s\n", A);
    printf("%s\n", B);

    // print_mat(F, m, n);

    int i = m - 2;
    int j = n - 2;

    // printf("\n%d_%d_%d_%d\n", m, n, i, j);

    char* AlignmentA = NULL;
    char* AlignmentB = NULL;
    AlignmentA = (char*)calloc(i, sizeof(char));
    AlignmentB = (char*)calloc(j, sizeof(char));

    while (i >= 0 && j >= 0) {
        if (i > 0 && j > 0 && F[i * m + j] == F[(i - 1) * n + j - 1] + (A[i] == B[j] ? match : mismatch)) {
            string_insert(AlignmentA, A + i);
            string_insert(AlignmentB, B + j);
            i--;
            j--;
        }
        else if (i > 0 && F[i * n + j] == F[(i - 1) * n + j] + gap) {
            string_insert(AlignmentA, A + i);
            string_insert(AlignmentB, "-");
            i--;
        }
        else {
            string_insert(AlignmentA, "-");
            string_insert(AlignmentB, B + j);
            j--;
        }
    }
    printf("%s\n", AlignmentA);
    printf("%s\n", AlignmentB);
}




int main(int argc, char** argv) {
    char A [] = "CGTGAATTCAT";
    char B [] = "GACTTAC";
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int match = 1;
    int mismatch = -1;
    int gap = -1;

    int* F = calculate_scoring_matrix(A, B, match, mismatch, gap);
    print_mat(F, m, n);

    align(F, A, B, match, mismatch, gap);
    // print_mat(F, m, n);

    free(F);
}